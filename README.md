# Projecting international mpox spread in Asia: ongoing global health risk

This repository contains the data and code for our paper:
> Toshiaki R. Asakura, Sung-mok Jung, Hiroaki Murayama, Cyrus Ghaznavi, Haruka
> Sakamoto, Ayaka Teshima, Fuminari Miura, Akira Endo (2024).
> Projecting international mpox spread in Asia:
> ongoing global health risk. medRxiv [Preprint].
> doi: [10.1101/2024.04.17.24305832](https://doi.org/10.1101/2024.04.17.24305832)

Note: Since the intermediate files for reproducing figures/tables are
very large and we cannot upload international flight volume data, some code will cause errors.
Therefore, this repository is prepared for the implementation part to be checked.

### Data used in our study
Our fitting process and simulations with meta-population model used the following data sets.
- Line list of mpox incidence in Japan: `/data/JPN_linelist/master_20230707.csv`.
- International flight volume data: `/data/flight/selected_flight_matrix.csv` (not uploaded).
- Population size in each country: `/data/pop_size_edit.csv`.

### How to run the code.
Clone this repository and type `docker compose up` to
install the Docker image and to set up the container.
Then, you can run the code in `script` via Jupyter Lab.

Run code in a sequence following the number of prefix of each file in `script` directory.
It takes around 1-2 days to complete each scenario in  `2j_AFP_MCMC.jl`, and this file will cause an error because of lacking international flight volume data.
`3j_prepare_data_for_vis.jl` prepared the necessary files for the visualisation. Outputs from `3j` file were uploaded in `tmp_results` and with thouse outputs, you can run `4j`, `5r`, `6r` and `7p` files to reproduce our results.

### License
[MIT](/LICENSE)
