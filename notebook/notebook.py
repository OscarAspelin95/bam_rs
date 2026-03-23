# Copyright 2026 Marimo. All rights reserved.

import marimo

__generated_with = "0.21.1"
app = marimo.App(css_file="style.css")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Picking File
    We'll use the built in file browser to access our .tsv file.
    """)
    return


@app.cell
def _(mo):
    file_picker = mo.ui.file_browser()

    file_picker
    return (file_picker,)


@app.cell(hide_code=True)
def _(mo, num_positions_input):
    mo.md(rf"""
    # Reading File
    We'll use polars to read a relatively large .tsv file. Specifically, we have aligned reads against a reference and have extracted information about coverage and deletion rates.

    To improve performance, only the first **{num_positions_input.value}** positions are used.
    """)
    return


@app.cell
def _(mo):
    num_positions_input = mo.ui.number(
        start=1,
        stop=1_000_000,
        step=1,
        value=50_000,
        debounce=True,
        label="Num positions to use",
    )

    num_positions_input
    return (num_positions_input,)


@app.cell
def _(file_picker, num_positions_input):
    from pathlib import Path

    import polars as pl

    def get_tsv() -> Path:
        match file_picker.value:
            case [file_info]:
                if not (tsv := file_info.path).is_file():
                    raise FileNotFoundError(f"{str(tsv)}")

                return file_info.path

            case _ as other:
                raise ValueError(f"Invalid input: `{other}`. Must be exactly one file.")

    _tsv = get_tsv()

    df = pl.read_csv(_tsv, separator="\t").head(num_positions_input.value)

    df
    return df, pl


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    # Plotting a DataFrame
    Here, we use altair to produce some nice, interactive plots.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualizing coverage
    We will plot sequencing depth as a function of position. By creating an interactive input, we can decide what data points we want to include. Change the **start**, **end** and **depth** thresholds to view the impact.
    """)
    return


@app.cell
def _(df, mo):
    _max_depth = df.max()["depth"][0]
    cov_input = mo.ui.dictionary(
        {
            "start": mo.ui.number(
                label="Start",
                value=1,
                start=0,
                step=1,
                stop=len(df),
                debounce=True,
            ),
            "end": mo.ui.number(
                label="End",
                value=min(1_000, len(df)),
                start=0,
                step=1,
                stop=len(df),
                debounce=True,
            ),
            "depth": mo.ui.range_slider(
                label="Depth",
                start=1,
                step=1,
                stop=_max_depth,
                show_value=True,
                debounce=True,
            ),
        }
    )
    return (cov_input,)


@app.cell
def _(cov_input):
    cov_input
    return


@app.cell
def _(cov_input, df, mo, pl):
    import altair as alt

    # alt.data_transformers.enable("vegafusion")

    _start = cov_input.elements["start"].value
    _end = cov_input.elements["end"].value
    _min_depth, _max_depth = cov_input.elements["depth"].value

    _df_alt = (
        df.filter(pl.col("pos").le(_end))
        .filter(pl.col("pos").ge(_start))
        .filter(pl.col("depth").ge(_min_depth))
        .filter(pl.col("depth").le(_max_depth))
    )

    selection = alt.selection_interval(encodings=["x"])

    tooltip = alt.Tooltip(field="depth")

    _chart = (
        alt.Chart(_df_alt, title="Sequencing Depth", width=500)
        .mark_line()
        .encode(x="pos", y="depth")
        .add_params(selection)
    )

    depth_chart = mo.ui.altair_chart(_chart)
    return alt, depth_chart


@app.cell
def _(depth_chart, mo):
    mo.vstack(
        [
            depth_chart,
            mo.md(
                "If the plot is empty, try bumping the **end** threshold to be larger than **start**."
            ),
        ]
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualizing base quality
    We have the mean base quality per reference base. We can visualize this to look for potential bias.
    """)
    return


@app.cell
def _(alt, df, mo):
    _chart = (
        alt.Chart(df, width=300, title="Base quality per reference base")
        .mark_boxplot()
        .encode(x="ref", y="base_phred", color="ref")
    )

    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Visualizing depth
    We have the depth per reference base. We can visualize this to look for potential bias.
    """)
    return


@app.cell
def _(alt, df, mo):
    _chart = (
        alt.Chart(df, width=300, title="Depth per reference base")
        .mark_boxplot()
        .encode(x="ref", y="depth", color="ref")
    )

    mo.ui.altair_chart(_chart)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Identifying Deletions
    Maybe, we find it interesting to find positions with a large deletion fraction. We can visualize this with a scatter plot of the deletion fraction against sequencing depth.
    """)
    return


@app.cell
def _(df, mo):
    _max_depth = df.max()["depth"][0]

    del_input = mo.ui.dictionary(
        {
            "min_frac_deletion": mo.ui.number(
                label="Minimum Deletion Fraction",
                value=0.5,
                start=0.1,
                step=0.1,
                stop=1.0,
                debounce=True,
            ),
            "depth": mo.ui.range_slider(
                label="Depth",
                start=1,
                step=1,
                stop=_max_depth,
                show_value=True,
                debounce=True,
            ),
        }
    )
    return (del_input,)


@app.cell
def _(del_input):
    del_input
    return


@app.cell
def _(alt, del_input, df, mo, pl):
    _min_frac_deletion = del_input.elements["min_frac_deletion"].value
    _min_depth, _max_depth = del_input.elements["depth"].value

    _alt_df = (
        df.filter(pl.col("frac_deletion").gt(_min_frac_deletion))
        .filter(pl.col("depth").ge(_min_depth))
        .filter(pl.col("depth").le(_max_depth))
    )

    _chart = (
        alt.Chart(_alt_df, width=550)
        .mark_point()
        .encode(x="depth", y="frac_deletion", color="ref")
    )

    mo.ui.altair_chart(_chart)
    return


if __name__ == "__main__":
    app.run()
