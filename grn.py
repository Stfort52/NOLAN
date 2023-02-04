import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests


def reconstructNet(
    net: pd.DataFrame,
    alpha: float = 0.05,
    apval: float | None = None,
    method: str = "fdr_bh",
    count: int | None = None,
    inplace: bool = False,
) -> pd.DataFrame:
    """
    Seralize and filter the TE matrix into a GRN.
    @params te: the TE matrix.
    @params alpha: the alpha value for false positive correction.
    @params apval: the adjusted p-value threshold.
    @params method: the method of false positive correction.
    @params count: the number of significant edges to be returned.
    """
    if not inplace:
        net = net.copy()
    relationship = net.to_numpy().flatten()
    source = net.index.repeat(len(net)).to_list()
    target = net.index.to_list() * len(net)

    regNet = pd.DataFrame(dict(source=source, relationship=relationship, target=target))
    regNet = regNet[regNet.source != regNet.target]
    mean, std = regNet.relationship.mean(), regNet.relationship.std()

    regNet["tstat"] = (regNet.relationship - mean) / std
    regNet["pval"] = 1 - stats.norm.cdf(regNet.tstat.astype(float))
    regNet["adj_pval"] = multipletests(regNet.pval, alpha=alpha, method=method)[1]

    regNet = regNet.sort_values("adj_pval")

    if apval is not None:
        regNet = regNet[regNet.adj_pval < apval]

    if count is not None:
        regNet = regNet.head(count)

    return regNet


def trimIndirect(
    net: pd.DataFrame, threshold: float = 0, depth: int = 1, inplace: bool = False
) -> pd.DataFrame:
    """
    Trim the indirect edges of a network.
    @params net: the filtered and Seralized network.
    @params threshold: the threshold of indirect edges.
    @params depth: the depth of indirect edges.
    """
    assert depth == 1, "Currently only supports depth=1"
    if not inplace:
        net = net.copy()
    net = net.sort_values("relationship")
    net['to_drop'] = False
    for idx, interaction in net.iterrows():
        ingress = net[net.source == interaction.source].drop(index=idx)
        egress = net[net.target == interaction.target].drop(index=idx)
        if not ingress.empty and not egress.empty:
            paths = pd.merge(
                ingress,
                egress,
                how="inner",
                left_on="target",
                right_on="source",
                suffixes=("_i", "_e"),
            )
            paths['relationship_min'] = paths[['relationship_i', 'relationship_e']].min(axis=1)
            if not paths[interaction.relationship < paths.relationship_min + threshold].empty:
                net.loc[idx, 'to_drop'] = True
    net = net[~net.to_drop]
    return net.drop(columns=['to_drop'])

