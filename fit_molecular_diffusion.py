import dataclasses
from typing import (
    Callable,
    Generator,
    Iterable,
    List,
    Tuple,
    TypeVar,
    Union,
    cast,
    overload,
)

import expression
import matplotlib.axes as mpl_axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns  # type: ignore
from statsmodels.regression.linear_model import OLS, RegressionResults  # type: ignore
from statsmodels.tools import add_constant  # type: ignore
from typing_extensions import Self

sns.set_style("whitegrid")  # type: ignore

a = TypeVar("a")
b = TypeVar("b")


def AICc(reg_res: RegressionResults) -> float:
    AIC: float = reg_res.aic  # type: ignore
    n: int = reg_res.nobs  # type: ignore
    k: int = reg_res.df_model + reg_res.k_constant  # type: ignore
    return AIC + 2 * k * (k + 1) / (n - k - 1)


@expression.curry(1)
def map(func: Callable[[a], b], iterable: Iterable[a]) -> Generator[b, None, None]:
    for x in iterable:
        yield func(x)


PFASs = pd.read_excel("PFASs.xlsx", sheet_name="PFASs").assign(
    n_C=lambda dataf: dataf["n_CFx"] + dataf["n_COO"] + dataf["n_CHx"],
    n_SO3xCFx=lambda dataf: dataf["n_SO3"] * dataf["n_CFx"],
)  # .query('PFAS != "FC-4"')  # type: ignore

x_vars = ["M", "n_C", "n_CFx", "n_COO", "n_SO3"]
y_vars = ["D0", "log D0"]

IndexType = TypeVar("IndexType", int, slice)
IndexType2 = TypeVar("IndexType2", int, slice)


class AxesArray1D:
    def __init__(self, axes: List[mpl_axes.Axes]) -> None:
        self.axes = axes

    @overload
    def __getitem__(self, key: int) -> mpl_axes.Axes:
        ...

    @overload
    def __getitem__(self, key: slice) -> "AxesArray1D":
        ...

    def __getitem__(self, key: IndexType) -> Union[mpl_axes.Axes, "AxesArray1D"]:
        if isinstance(key, slice):
            return AxesArray1D(self.axes[key])
        return self.axes[key]

    def __iter__(self) -> Generator[mpl_axes.Axes, None, None]:
        yield from self.axes


class AxesArray2D:
    def __init__(self, axes: List[List[mpl_axes.Axes]]) -> None:
        self.axes = axes

    @overload
    def __getitem__(self, key: Tuple[int, int]) -> mpl_axes.Axes:
        ...

    @overload
    def __getitem__(self, key: int) -> AxesArray1D:
        ...

    @overload
    def __getitem__(self, key: Tuple[slice, int]) -> AxesArray1D:
        ...

    @overload
    def __getitem__(self, key: Tuple[int, slice]) -> AxesArray1D:
        ...

    @overload
    def __getitem__(self, key: slice) -> "AxesArray2D":
        ...

    @overload
    def __getitem__(self, key: Tuple[slice, slice]) -> "AxesArray2D":
        ...

    def __getitem__(
        self, key: Union[IndexType, Tuple[IndexType, IndexType2]]
    ) -> Union[mpl_axes.Axes, AxesArray1D, "AxesArray2D"]:
        if isinstance(key, int):
            return AxesArray1D(self.axes[key])
        if isinstance(key, slice):
            return AxesArray2D(self.axes[key])
        if len(key) == 2:
            row, col = key
            if isinstance(row, int) and isinstance(col, int):
                return self.axes[row][col]
            if isinstance(row, slice) and isinstance(col, int):
                return AxesArray1D([ax[col] for ax in self.axes[row]])
            if isinstance(row, int) and isinstance(col, slice):
                return AxesArray1D(self.axes[row][col])
            if isinstance(row, slice) and isinstance(col, slice):
                return AxesArray2D([ax[col] for ax in self.axes[row]])
        raise NotImplementedError

    def __iter__(self) -> Generator[AxesArray1D, None, None]:
        for row in self.axes:
            yield AxesArray1D(row)


fig, axs_ = plt.subplots(  # type: ignore
    nrows=len(y_vars),
    ncols=len(x_vars),
    sharex="col",
    sharey="row",
    figsize=(1 + 2.5 * len(x_vars), 1 + 2.5 * len(y_vars)),
    dpi=200,
    constrained_layout=True,
)
axs = AxesArray2D(axs_)

for i, x_var in enumerate(x_vars):
    for j, y_var in enumerate(y_vars):
        sns.regplot(  # type: ignore
            data=PFASs,
            x=x_var,
            y=y_var,
            ax=axs[j, i],
            scatter_kws={"s": 10, "alpha": 0.5},
            logx=False,
        )

for label, ax in zip(y_vars, axs[:, 0]):
    ax.set_ylabel(label, fontsize=12)  # type: ignore
for label, ax in zip(x_vars, axs[-1, :]):
    ax.set_xlabel(label, fontsize=12)  # type: ignore
fig.suptitle("Linear regression of molecular diffusion on molecular properties")  # type: ignore
fig.savefig("molecular_diffusion.png", dpi=200)  # type: ignore
plt.close(fig)  # type: ignore


@dataclasses.dataclass
class ModelVariant:
    name: str = dataclasses.field(repr=True)
    variables: pd.DataFrame = dataclasses.field(repr=False)


@dataclasses.dataclass
class FittedModelVariant(ModelVariant):
    model: RegressionResults = dataclasses.field(repr=False)
    AICc: float = dataclasses.field(repr=True)


@dataclasses.dataclass
class ModelSet:
    variants: List[FittedModelVariant]

    def __post_init__(self) -> None:
        self.variants.sort(key=lambda v: v.AICc)
        self.best_model = self.variants[0]
        self.relative_performance = [
            np.exp(0.5 * (self.best_model.AICc - v.AICc)) for v in self.variants
        ]

    @classmethod
    def from_variants(cls, variants: Iterable[FittedModelVariant]) -> Self:
        return ModelSet(list(variants))

    def __len__(self) -> int:
        return len(self.variants)


y: "pd.Series[float]" = PFASs["log D0"]

model_variants = [
    ModelVariant("M", PFASs[["M"]]),
    ModelVariant("#CFx", PFASs[["n_CFx"]]),
    ModelVariant("#CFx and #SO3*#CFx", PFASs[["n_CFx", "n_SO3xCFx"]]),
    ModelVariant("M and #CFx", PFASs[["M", "n_CFx"]]),
    ModelVariant("M, #CFx and #COO", PFASs[["M", "n_CFx", "n_COO"]]),
    ModelVariant(
        "M, #CFx, #COO and #SO3",
        PFASs[["M", "n_CFx", "n_COO", "n_SO3"]],
    ),
]


@expression.curry_flipped(1)
def fit_model_variant(
    model_variant: ModelVariant, model_target: "pd.Series[float]"
) -> FittedModelVariant:
    X: pd.DataFrame = add_constant(model_variant.variables)  # type: ignore
    model = cast(RegressionResults, OLS(model_target, X).fit())  # type: ignore
    return FittedModelVariant(
        model_variant.name, model_variant.variables, model, AICc(model)
    )


model_set = expression.pipe(
    model_variants,
    map(fit_model_variant(y)),
    ModelSet.from_variants,
)
print(model_set.best_model.model.summary())
print("Selected model:", model_set.best_model.name)
print("Intercept:", model_set.best_model.model.params["const"])
print("Slope:", model_set.best_model.model.params["n_CFx"])


fig, axs = plt.subplots(
    ncols=len(model_set),
    sharey=True,
    sharex=True,
    figsize=(1 + 2.5 * len(model_set), 3.5),
    dpi=200,
    constrained_layout=True,
)
for ax, model_variant in zip(axs, model_set.variants):
    sns.scatterplot(
        x=y.values,
        y=model_variant.model.fittedvalues,
        hue=PFASs["n_CFx"].values,
        style=PFASs["n_SO3xCFx"].values,
        ax=ax,
    )
    sns.regplot(x=y.values, y=model_variant.model.fittedvalues, ax=ax, scatter=False)
    ax.axline(xy1=(y.mean(), y.mean()), slope=1, color="k", linestyle="--")
    ax.set_aspect("equal")
    ax.set_title(model_variant.name)
ax.set_xlim(ax.get_ylim())
fig.supxlabel("log D0")
fig.supylabel("log D0 (predicted)")
fig.suptitle("Linear regression of molecular diffusion on molecular properties")
fig.savefig("molecular_diffusion_linear_regression.png", dpi=200)
