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

# # OWID data up to the end of 2024

import datetime
import polars as pl
import matplotlib.pyplot as plt

df = pl.read_csv("../data/owid-mpox-data-20250121.csv")
df = df.with_columns(
    pl.col("date").str.to_datetime(format="%d/%m/%Y")
)

# # Asia 2022

df2022 = df.filter(df["date"] < datetime(2023, 1, 1))
df_global = df2022.filter(df2022["location"] == "World")
df_Asia = df2022.filter(df2022["location"] == "Asia")
df_Israel = df2022.filter(df2022["location"] == "Israel")

n2022_Asia = df_Asia.tail(1)[0, "total_cases"]
n2022_Israel = df_Israel.tail(1)[0, "total_cases"]
print(f"{n2022_Israel / n2022_Asia * 100:.2f}% ({n2022_Israel}/{n2022_Asia})")

plt.plot(df_Asia["date"], df_Asia["total_cases"])
plt.plot(df_Israel["date"], df_Israel["total_cases"])
plt.xlim([datetime(2022, 5, 1), datetime(2022,12,31)])

# # Cumulative incidence for each year

Asia_list = pl.read_csv("../data/pop_size_edit.csv")
iso_cnts = Asia_list["iso_code"].to_list()

dfM = df.filter( df["iso_code"].is_in(iso_cnts) )
cond1 = (pl.col("new_cases") > 0 ) & \
        (pl.col("date") >= datetime.datetime(2022, 1, 1)) & \
        (pl.col("date") <= datetime.datetime(2022, 12,31))
cond2 = (pl.col("new_cases") > 0 ) & \
        (pl.col("date") >= datetime.datetime(2023, 1, 1)) & \
        (pl.col("date") <= datetime.datetime(2023, 12,31))
cond3 = (pl.col("new_cases") > 0 ) & \
        (pl.col("date") >= datetime.datetime(2024, 1, 1)) & \
        (pl.col("date") <= datetime.datetime(2024, 12,31))


res = dfM.group_by("location").agg(
                                  pl.col("new_cases").filter(cond1)
                                      .sum().alias("cumsum_2022"),
                                  pl.col("new_cases").filter(cond2)
                                      .sum().alias("cumsum_2023"),
                                  pl.col("new_cases").filter(cond3)
                                      .sum().alias("cumsum_2024"),
                                  pl.col("date").filter(cond2)
                                      .min().alias("first_report_2023"),
                                  pl.col("date").filter(cond2)
                                      .max().alias("latest_date_2023"),
                                 )
res = res.sort(by="first_report_2023", descending=False)
res = res.with_columns(
    pl.col("latest_date_2023").dt.strftime(format="%Y-%m-%d"),
    pl.col("first_report_2023").dt.strftime(format="%Y-%m-%d")
)
#res = Asia_list[["location"]].join(res, on="location", how="outer")
res.write_csv("../fig/owid_summary_table.csv")
res

lis = ["cumsum_2023", "cumsum_2024"]
res[lis].sum()


