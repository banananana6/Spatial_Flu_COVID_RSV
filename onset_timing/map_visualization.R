# visualization of breakpoints on U.S. map

library(ggplot2)
library(dplyr)
library(maps)
library(mgcv)

# start_date <-"2023-03-20"
# end_date <- "2026-03-01"

us_map <- map_data("state")

hsa_geo <- read.csv("HSA_GEO.csv")
breakpoints <- read.csv("breakpoints_covid1.csv")
breakpoints <- breakpoints %>%
  mutate(breakpoint = as.Date(breakpoint, format = "%m/%d/%Y")) %>%
  # mutate(breakpoint = as.Date(breakpoint, format = "%Y-%m-%d")) %>%
  # filter(breakpoint >= as.Date(start_date) & breakpoint <= as.Date(end_date)) %>%
  filter(hsa_nci_id != "All") %>%
  group_by(hsa_nci_id) %>%
  mutate(hsa_nci_id=as.numeric(hsa_nci_id)) %>%
  # slice_min(breakpoint, n = 1, with_ties = FALSE) %>%
  ungroup()

## --- to filter out outliers ---

# bp_numeric <- as.numeric(breakpoints$breakpoint)
# 
# q1 <- quantile(bp_numeric, 0.25, na.rm = TRUE)
# q3 <- quantile(bp_numeric, 0.75, na.rm = TRUE)
# iqr <- q3-q1

# breakpoints <- breakpoints |>
#   filter(
#     as.numeric(breakpoint) >= q1-1.5*iqr,
#     as.numeric(breakpoint) <= q3+1.5*iqr
#   )

hsa_bp <- hsa_geo %>%
  left_join(breakpoints, by = c("hsa_id" = "hsa_nci_id")) %>%
  filter(!is.na(breakpoint))
hsa_bp <- hsa_bp %>% 
  mutate(breakpoint = as.Date(breakpoint, format = "%m-%d-%Y"))

bp_counts <- breakpoints %>%
  group_by(hsa_nci_id) %>%
  summarise(n_breakpoints = n()) %>%
  left_join(hsa_geo, by = c("hsa_nci_id" = "hsa_id")) %>%
  filter(!is.na(weighted_lon), !is.na(weighted_lat))

# --- point visualization on map ---

# ggplot() +
#   coord_map("albers", lat0 = 39, lat1 = 45) +
#   geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
#                fill = "grey90", color = "white", linewidth = 0.3) +
#   geom_point(data = hsa_bp, aes(x = weighted_lon, y = weighted_lat, color = breakpoint),
#              size = 2, alpha = 0.8) +
#   scale_color_date(name = "Breakpoint", date_labels = "%b %Y", low = "red", high = "blue")+
#   theme_void(base_family = "") +
#   theme(legend.position = "right")

## --- number of breakpoints, plotted on map ---

# ggplot() +
#   coord_map("albers", lat0 = 39, lat1 = 45) +
#   geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
#                fill = "grey90", color = "white", linewidth = 0.3) +
#   geom_point(data = bp_counts, aes(x = weighted_lon, y = weighted_lat, color = n_breakpoints),
#              size = 2, alpha = 0.8) +
#   scale_color_viridis_c(option = "plasma", name = "# Breakpoints") +
#   theme_void(base_family = "") +
#   theme(legend.position = "right")

# --- smooth interpolation of breakpoints ---

grid <- expand.grid(
  weighted_lon = seq(-125, -65, length.out = 200),
  weighted_lat = seq(24, 50, length.out = 200)
)

hsa_bp$breakpoint_num <- as.numeric(hsa_bp$breakpoint)
fit <- gam(breakpoint_num ~ s(weighted_lon, weighted_lat), data = hsa_bp)
grid$pred_num <- predict(fit, newdata = grid)
grid$pred_date <- as.Date(grid$pred_num, origin = "1970-01-01")

ggplot() +
  geom_tile(data = grid, aes(x = weighted_lon, y = weighted_lat, fill = pred_date)) +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_date(low = "blue", high = "red", name = "Breakpoint Date") +
  coord_fixed(1.3, xlim = c(-125, -65), ylim = c(24, 50)) +
  theme_void() +
  labs(title = "Smoothed COVID Breakpoint Dates across U.S. HSAs")

# --- smooth interpolation of breakpoint frequency ---

grid <- expand.grid(
  weighted_lon = seq(-125, -65, length.out = 200),
  weighted_lat = seq(24, 50, length.out = 200)
)

fit <- gam(n_breakpoints ~ s(weighted_lon, weighted_lat), data = bp_counts)
grid$pred_num <- predict(fit, newdata = grid)

ggplot() +
  geom_tile(data = grid, aes(x = weighted_lon, y = weighted_lat, fill = pred_num)) +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = NA, color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "# Breakpoints", option = "plasma")+
  coord_fixed(1.3, xlim = c(-125, -65), ylim = c(24, 50)) +
  theme_void() +
  labs(title = "Smoothed RSV Breakpoint Frequency across U.S. HSAs")
