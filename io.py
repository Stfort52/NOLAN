import pandas as pd
import numpy as np
import os


def parseGenes(
    path: str | os.PathLike, genesOnRows: bool = False, delim: str = ","
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """
    parse Wishbone flavored single cell expression matrix
    @params path: the path to the file.
    @params genesOnRows: if the genes are on the row indices.
    @params delim: the delimiter of the file.
    @return a tuple of (expression, cell label, genes)
    """
    df: pd.DataFrame = pd.read_csv(path, delimiter=delim)
    df.rename_axis(columns=None, index=None, inplace=True)
    if genesOnRows:
        df = df.T
        df, df.columns = df.iloc[1:], df.iloc[0]  # row 0 as columns
    else:
        df, df.index = df.iloc[:, 1:], df.iloc[:, 0]
    cells = df.index.to_series()
    genes = df.columns.to_series()
    return df, cells, genes


def parseCellSelect(path: str | os.PathLike) -> list[int]:
    with open(path) as f:
        return list(map(int, f))


def parseTrajectory(path: str | os.PathLike) -> list[float]:
    with open(path) as f:
        return list(map(float, f))


def writeTEmtx(path: str | os.PathLike, te: pd.DataFrame):
    with open(path, "w") as f:
        f.write("TE")
        te.to_csv(f, sep="\t")


def readTEmtx(path: str | os.PathLike) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0).rename_axis(index=None)

def locateJar(path: str|os.PathLike|None = None) -> str:
    """
    locate the jar file
    @params path: the path to the jar file or the directory containing the jar file.
    @return the path to the jar file.
    Default path is the submodule's jar.
    """
    if not path:
        return os.path.join(os.path.dirname(__file__), "TENET",  "infodynamics.jar")
    elif os.path.basename(path) == "infodynamics.jar" and os.path.exists(path):
        return path
    elif os.path.isdir(path) and os.path.exists(os.path.join(path, "infodynamics.jar")):
        return os.path.join(path, "infodynamics.jar")
    else:
        raise FileNotFoundError(f"{path} not found")
