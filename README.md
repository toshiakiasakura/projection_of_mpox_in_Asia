# Projecting international mpox spread in Asia: ongoing global health risk

This repository contains the data and code for our study.
<!--- for our paper: [Add citation info] -->

Note: Since the intermediate files for reproducing figures/tables are
very large and we cannot upload international flight volume data, some code will cause errors.
Therefore, this repository is prepared to share the implementation code.

### Data used in our study
Our fitting process and simulations with meta-population model used the following data sets.
- Line list of mpox incidence in Japan: `/data/JPN_linelist/master_20230707.csv`.
- International flight volume data: `/data/flight/selected_flight_matrix.csv` (not uploaded).
- Population size in each country: `/data/pop_size_edit.csv`.

### How to run the code.
Clone this repository and type `docker compose up` to
install the Docker image and to set up the container.
Then, you can run the code via Jupyter Lab.

The files were executed in a sequence following the number of prefix of each file name in `script` directory.
It takes around 1-2 days to complete each scenario in  `2j_AFP_MCMC.ipynb` and this code can not be executed because of the lack of international flight volume data.
`3j_prepare_data_for_vis.jl` prepared the necessary files for the visualisation. Outputs from this file were uploaded in thie repository and with thouse outputs, you can run `4j`, `5r`, `6r` and `7p` to reproduce our results.

