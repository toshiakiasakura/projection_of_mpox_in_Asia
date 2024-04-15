# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import datetime
import polars as pl

df = pl.read_csv("../data/owid-mpox-data-20240303.csv")
df = df.with_columns(
    pl.col("date").str.to_datetime(format="%Y-%m-%d")
)

Asia_list = pl.read_csv("../data/pop_size_edit.csv")

iso_cnts = Asia_list["iso_code"].to_list()

dfM = df.filter( df["iso_code"].is_in(iso_cnts) )
cond1 = (pl.col("new_cases") > 0 ) & \
        (pl.col("date") >= datetime.datetime(2023, 1, 1)) & \
        (pl.col("date") <= datetime.datetime(2023, 12,31))
cond2 = (pl.col("new_cases") > 0 ) & \
        (pl.col("date") < datetime.datetime(2023, 1, 1))

res = dfM.group_by("location").agg(
                                  pl.col("new_cases").filter(cond2)
                                      .sum().alias("cumsum_2022"),
                                  pl.col("new_cases").filter(cond1)
                                      .sum().alias("cumsum_2023"),
                                  pl.col("date").filter(cond1)
                                      .max().alias("latest_date_2023"),
                                  pl.col("date").filter(cond1)
                                      .min().alias("first_report_2023")
                                 )
res = res.sort(by="first_report_2023", descending=False)
res = res.with_columns(
    pl.col("latest_date_2023").dt.strftime(format="%Y-%m-%d"),
    pl.col("first_report_2023").dt.strftime(format="%Y-%m-%d")
)
#res = Asia_list[["location"]].join(res, on="location", how="outer")
res.write_csv("../fig/owid_summary_table.csv")
res





