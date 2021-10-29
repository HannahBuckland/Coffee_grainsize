# Script to process raw CX2 files
# Hannah Buckland 31/07/2021

# Loading required R packages

library(tidyverse)
library(data.table)
library(patchwork)
library(viridis)

######################################################################
# First - Read in .xle files (excel compatible files that the Camsizer produces)


# List all files witl .xle extension
all_xle <-
  list.files(
    "Camsizer_raw",
    pattern = ".xle",
    recursive = TRUE,
    full.names = TRUE
  )

# Extract ids for each analysis from file name
samplename <-
  str_extract(all_xle, "(?<=Camsizer_raw/).*(?=/)") # isolates the sample name
analysisid <- str_extract(all_xle, "(?<=Camsizer_raw/).*(?=.xle)")
analysisid <-
  tstrsplit(analysisid, "/")[[2]] # then we get the long analysis name (there are repeat runs for multiple samples and different size parameters)

# need to isolate the size parameter for each analysis (see Buckland et al. 2021 for definitions)
paramid <- ifelse(grepl("x_area", analysisid),
                  "x_area",
                  ifelse(grepl("xc_min", analysisid), "xc_min",
                         "xFe_max"))

# designating simplified name for each analysis (will have replicates for repeated analyses)
simple_name <-
  paste0(samplename, "_", paramid)


# Function to read in .xle files
file_process <- function(xlefile) {
  # need to skip nulls and set the file encoding because of strange format output directly from the CX2
  dat <- read.table(
    xlefile,
    fill = TRUE,
    skipNul = TRUE,
    fileEncoding = "latin1",
    stringsAsFactors = FALSE
  )
  
  # only extracting the distributions
  dat_clean <- tail(dat, 98)
  # ignoring suprious coloum
  dat_clean <- dat_clean[, 1:8]
  # set column names
  colnames(dat_clean) <-
    c("min_um",
      "max_um",
      "vol_perc",
      "cum1",
      "SPHT",
      "Symm",
      "bl",
      "PND")
  
  # here need to convert data into numerics
  dat_clean <- dat_clean %>%
    mutate_if(is.character, as.numeric) %>%
    arrange(desc(min_um)) %>% # sort by descending size
    mutate(cum2 = cumsum(vol_perc)) %>% # recalcualte cumulative sum as coarser than percentage
    arrange(min_um) %>% # put back to original order
    mutate(min_um = ifelse(min_um == 0, 0.8,
                           min_um)) # set minimum particle size of CX2 instrument which is 0.8 µm (see Buckland et al. 2021)
  
  return(dat_clean) # return the cleaned data from the function
  
}

# Apply the function across all the .xle files to read them in and clean them up
all_xlefiles <- lapply(all_xle, FUN = file_process)
# Assign names to the all_xlefiles list and unique=TRUE allows for repeat measurements to be assigned a number
names(all_xlefiles) <- make.names(simple_name,
                                  unique = TRUE)

# Unlist the list into long format where the id column is equal to the name of the list item
xlefiles_bind <- rbindlist(all_xlefiles, idcol = TRUE)


# Now need to add a couple of parameters to each run (sample id, size parameter, bean type and grind setting)
xlefiles_processed <- xlefiles_bind %>%
  mutate(
    sample = tstrsplit(.id, "_")[[1]],
    param = paste0(tstrsplit(.id, "_")[[2]], tstrsplit(.id, "_")[[3]]),
    param = gsub("([[:punct:]]+)([[:digit:]]+)", '', param),
    name = paste0(sample, param),
    bean = ifelse(grepl("Guayancan", sample), "Guayancan",
                  "Tumba"),
    grind = paste0("grind", tstrsplit(sample, bean)[[2]]),
    #want to change grind variable to an ordered factor so that I can assign colours in plotting
    grind = factor(
      grind,
      levels = c(
        "grind01",
        "grind02",
        "grind03",
        "grind04",
        "grind05",
        "grind06",
        "grind07",
        "grind08",
        "grind09",
        "grind10",
        "grind11"
      )
    )
  ) %>%
  select(-.id)

######################################################################
# Dealing with repeat analyses

# Find out where the duplicated analyses are
xle_duplicates <- xlefiles_processed %>%
  group_by(name) %>%
  summarise(count = n())

# Get a vector of each run name that we want an average for
unique_runs <- unique(xle_duplicates$name)

# Function to average the repeat runs for each sample and size parameter combo
average_repeats <- function(runname) {
  subset_xle <- xlefiles_processed %>%
    filter(name == runname) %>%
    select(c(
      name,
      sample,
      param,
      grind,
      bean,
      min_um,
      max_um,
      vol_perc,
      cum2,
      SPHT,
      bl,
      PND
    )) %>%
    pivot_wider(
      names_from =  name,
      values_from = c(vol_perc, cum2, SPHT, bl, PND),
      values_fn = mean
    )
  
  # can't figure out how to now have the prefix on the variables so manually renaming columns
  colnames(subset_xle) <-
    c(
      "sample",
      "param",
      "grind",
      "bean",
      "min_um",
      "max_um",
      "vol_perc",
      "cum2",
      "SPHT",
      "bl",
      "PND"
    )
  return(subset_xle)
  
}
# Apply the averaging function
averaged_xle <- lapply(unique_runs, FUN = average_repeats)
# Assign names to the output of the function
names(averaged_xle) <- make.names(unique_runs)

#Unlist the list into long format
averaged_bind <- rbindlist(averaged_xle, idcol = TRUE)

# Final check will show that there is now only one entry dataset (an average) for each grind setting and size parameter
xle_no_dups <- averaged_bind %>%
  group_by(.id) %>%
  summarise(count = n())

######################################################################
# Plotting code below this line

# Setting up colour palette in coffee(ish) colours
pallete <- c(
  grind01 = "#edc4b3",
  grind02 = "#e6b8a2",
  grind03 = "#deab90",
  grind04 = "#d69f7e",
  grind05 = "#cd9777",
  grind06 = "#c38e70",
  grind07 = "#b07d62",
  grind08 = "#9d6b53",
  grind09 = "#8a5a44",
  grind10 = "#774936",
  grind11 = "#1d0200"
)


# Produce plot to demonstrate how grind setting relates to the actual GSD

Cplots <- ggplot(data = averaged_bind %>%
                   filter(param %in% c("xFemax", "xcmin"),
                          bean == "Tumba")) +
  geom_hline(yintercept = 50,
             colour = "grey80",
             linetype = "dotted") +
  geom_line(aes(
    x = min_um,
    y = cum2,
    colour = grind,
    linetype = bean
  )) +
  scale_colour_manual(values = pallete, labels = seq(1, 11, by = 1)) +
  scale_linetype_manual(values = c(1)) +
  scale_x_log10(name = "Grain size (µm)") +
  scale_y_continuous(name = "Cumulative Volume %", expand = c(0.01, 0.01)) +
  labs(colour = "Grind", linetype = "Bean type") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    axis.text = element_text(colour = "black", size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Facet by the size parameter and give title
grindsize <-
  Cplots + facet_wrap( ~ param) + ggtitle("Grind size to grain size")

# Produce plot to demonstrate how the bean type influences the GSD produced by the same grind setting

Bplots <- ggplot(data = averaged_bind %>%
                   filter(param %in% c("xFemax", "xcmin"))) +
  geom_hline(yintercept = 50,
             colour = "grey80",
             linetype = "dotted") +
  geom_line(aes(
    x = min_um,
    y = cum2,
    colour = grind,
    linetype = bean
  )) +
  scale_colour_manual(values = pallete, labels = seq(1, 11, by = 1)) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_x_log10(name = "Grain size (µm)") +
  scale_y_continuous(name = "Cumulative Volume %", expand = c(0.01, 0.01)) +
  labs(colour = "Grind", linetype = "Bean type") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    axis.text = element_text(colour = "black", size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# Facet by the size parameter and give title
beantype <-
  Bplots + facet_wrap( ~ param) + ggtitle("Comparing bean types")

# Save individual plots as pngs
ggsave(
  "plots/grind_size.png",
  grindsize,
  width = 8,
  height = 4,
  units = "in"
)
ggsave(
  "plots/beantype.png",
  beantype,
  width = 8,
  height = 4,
  units = "in"
)
