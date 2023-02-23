from NOLAN.parallel import SharedArray, roundDown2sqrt, splitRanges
from NOLAN.io import locateJar
from itertools import combinations, product
import pandas as pd
import multiprocessing as mp
from jpype import *
from tqdm import tqdm
import numpy as np
import os

def calcTE(
    expr: pd.DataFrame,
    xRange: tuple,
    yRange: tuple,
    result: np.ndarray | SharedArray,
    histLen: int = 1,
    kernel: str|os.PathLike = locateJar(),
):
    """
    Calculates the transfer entropy value over the given range on the expression matrix.
    This function is for internal use, and not imported on the top level.
    """
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=" + kernel, "-Xmx16G")
    kern = JPackage(
        "infodynamics.measures.continuous.kernel"
    ).TransferEntropyCalculatorKernel
    if isinstance(result, SharedArray):
        result = result.to_numpy()

    for i in range(*xRange):
        for j in range(*yRange):
            if i == j:
                result[i, j] = 0
            else:
                teCalc = kern()
                teCalc.setProperty("NORMALISE", "true")
                teCalc.initialise(histLen, 0.5)
                teCalc.setObservations(
                    JArray(JDouble, 1)(expr.iloc[:, i].values),
                    JArray(JDouble, 1)(expr.iloc[:, j].values),
                )
                te = teCalc.computeAverageLocalOfObservations()
                result[i, j] = te


def inferNetwork(
    expr: pd.DataFrame,
    genes: list[str] | pd.Series,
    trajectory: list[float] | pd.Series,
    mask: list[bool] | pd.Series,
    histLen: int = 1,
    workers: int | tuple[int] = 1,
    kernel: str|os.PathLike = locateJar(),
) -> pd.DataFrame:
    """
    Infer the network from the expression data and the trajectory.
    @params `expr` expression data, a pandas.DataFrame. Genes on Columns, Cells on Rows.
    @params `genes` the subset of genes to use in the analysis.
    @params `trajectory` the trajectory, a Series of cell pseudotime.
    @params `mask` the mask of cells to be used, a Series of bool.
    @params `histLen` the history length to use in the TE calculation step.
    @params `workers` the number of workers to use, or tuple of two integers for block height and width.
    @params `kernel` the java archive containing the transfer entropy calcuation kernel.

    NOLAN will use the original TENET's `information.jar` as the default kernel in the TENET submodule.
    The original TENET paper recommends the history length of 1 for the best results.
    """

    def runTE(
        x1: np.ndarray,
        x2: np.ndarray,
        histLen: int,
        kern: JClass
    ):
        teCalc = kern()
        teCalc.setProperty("NORMALISE", "true")
        teCalc.initialise(histLen, 0.5)
        teCalc.setObservations(JArray(JDouble, 1)(x1), JArray(JDouble, 1)(x2))
        return teCalc.computeAverageLocalOfObservations()

    if isinstance(trajectory, list):
        trajectory = pd.Series(trajectory)

    if isinstance(mask, list):
        mask = pd.Series(mask)

    assert len(trajectory) == len(mask), "Length of trajectory and mask must be equal"
    assert len(trajectory) == len(
        expr.index
    ), "Length of trajectory and number of cells must be equal"

    # mask and sort expression matrix
    expr = expr.assign(
        trajectory=trajectory.values, cellSelect=mask.astype(bool).values
    )
    expr = expr[expr.cellSelect]
    expr = expr.sort_values(by="trajectory")
    expr = expr.drop(columns=["trajectory", "cellSelect"])

    ## Sequential Code
    if isinstance(workers, int) and workers <= 1:
        if not isJVMStarted():
            startJVM(
                getDefaultJVMPath(),
                "-ea",
                "-Djava.class.path=" + kernel,
                "-Xmx16G",
            )

        kern = JPackage(
            "infodynamics.measures.continuous.kernel"
        ).TransferEntropyCalculatorKernel
        result = pd.DataFrame(index=genes, columns=genes)

        for i in range(len(genes)):
            result.iloc[i, i] = 0

        for comb in tqdm(combinations(range(len(genes)), 2)):
            x1 = expr.iloc[:, comb[0]].to_list()
            x2 = expr.iloc[:, comb[1]].to_list()
            r1, r2 = runTE(x1, x2, histLen, kern), runTE(x2, x1, histLen, kern)
            result.iloc[comb[0], comb[1]] = r1
            result.iloc[comb[1], comb[0]] = r2

    else:
        assert (
            not isJVMStarted()
        ), "Forking Processes after the start of Jpype JVM breaks the JVM"
        # About: https://github.com/jpype-project/jpype/issues/1024

        if isinstance(workers, tuple):
            assert len(workers) == 2, "workers must be a tuple of length 2"
            blocks = product(
                splitRanges(len(genes), workers[0]), splitRanges(len(genes), workers[1])
            )
        else:
            workers = roundDown2sqrt(workers)
            blocks = product(splitRanges(len(genes), workers), repeat=2)
        shm = SharedArray(float, (len(genes), len(genes)))

        procs:list[mp.Process] = []

        for i, blk in enumerate(blocks):
            proc = mp.Process(
                target=calcTE,
                args=(expr, *blk, shm, histLen, kernel),
                name=f"calcTE_{i}@({blk})",
            )
            proc.start()
            procs.append(proc)

        for proc in procs:
            proc.join()

        result = pd.DataFrame(shm.to_numpy(), index=genes, columns=genes)

    return result
