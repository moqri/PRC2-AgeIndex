library(ClockBasis)
library(methylCIPHER)
tidylog()

source("/Lasagna/kying/p_2022/GEO_data/src/4a_compute_causal_age.R")
source("/Lasagna/kying/p_2022/GEO_data/src/4b_compute_grimage.R")

## Read data

files = dir('./data/', pattern = '*.csv')
data = read_csv(paste0('./data/', files), id = "file")

fdata = data |>
  mutate(file = basename(file) |> str_remove(".csv")) |>
  select(file, Probe_ID, beta) |>
  pivot_wider(names_from = file, values_from = beta)

normal_clocks = do_dnam_clock_human(fdata)

causal_clocks = causal_clock(fdata)

grim_age = do_grimage1(fdata, keep = F)

grim_age2 = do_grimage2(fdata, keep = T)

result = bind_cols(normal_clocks, causal_clocks[, -1], grim_age[, -1], grim_age2[, -1]) |>
  select(-Age, -Female, -ends_with('missing'))

write_csv(result, "./all_clock_res.csv")


## Read new data
files = dir('./data/s_arrays', pattern = '*.csv')
data = read_csv(paste0('./data/s_arrays/', files), id = "file")

fdata = data |>
  mutate(file = basename(file) |> str_remove(".csv")) |>
  select(file, Probe_ID, beta) |>
  pivot_wider(names_from = file, values_from = beta)

normal_clocks = do_dnam_clock_human(fdata)

causal_clocks = causal_clock(fdata)

grim_age = do_grimage1(fdata, keep = F)

grim_age2 = do_grimage2(fdata, keep = T)

result = bind_cols(normal_clocks, causal_clocks[, -1], grim_age[, -1], grim_age2[, -1]) |>
  select(-Age, -Female, -ends_with('missing'))

write_csv(result, "./all_clock_res2.csv")

## Read new data
files = dir('./data/f_arrays', pattern = '*.csv')
data = read_csv(paste0('./data/f_arrays/', files), id = "file")

fdata = data |>
  mutate(file = basename(file) |> str_remove(".csv")) |>
  select(file, Probe_ID, beta) |>
  pivot_wider(names_from = file, values_from = beta)

normal_clocks = do_dnam_clock_human(fdata)

causal_clocks = causal_clock(fdata)

grim_age = do_grimage1(fdata, keep = F)

grim_age2 = do_grimage2(fdata, keep = T)

result = bind_cols(normal_clocks, causal_clocks[, -1], grim_age[, -1], grim_age2[, -1]) |>
  select(-Age, -Female, -ends_with('missing'))

write_csv(result, "./f_clocks.csv")
