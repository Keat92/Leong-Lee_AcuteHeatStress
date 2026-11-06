# =============================================================================
# Part 1 — 2024 marine heatwave: in-situ (MESN buoy) and satellite SST / DHW
# =============================================================================
# Quantifies thermal stress during the 2024 marine heatwave from:
#   * Temperature data from the MESN buoy
#     (March 2023 – September 2024), and
#   * NOAA Coral Reef Watch CoralTemp v3.1 daily 5 km SST and Degree Heating
#     Weeks (DHW) for the pixel nearest the buoy.
#
# Site Maximum Monthly Mean (MMM) = 29.9 °C (supplied site value).
# Bleaching threshold = MMM + 1 °C (NOAA CRW HotSpot definition).
#
# Inputs  (data/raw/):
#   mesn_data.csv             buoy record: timestamp, temp, qcflag
#   noaa_coraltemp_2024.csv   CoralTemp SST + DHW, 2024 (downloaded if absent)
#
# =============================================================================

library(tidyverse)
library(lubridate)
library(zoo)        # rollapplyr() for the 84-day DHW window
library(scales)
library(rerddap)    # NOAA ERDDAP download (only used if data file is absent)
library(patchwork)  # inset_element()
library(ggpubr)     # theme_pubr()
library(here)

# -----------------------------------------------------------------------------
# 1. Load and clean buoy data
# -----------------------------------------------------------------------------
# The raw file mixes 1-minute and 5-minute samples, contains NA rows between
# buoy deployments, and carries a QARTOD-style qcflag (1 = good, 3 = suspect).
raw <- read_csv(here("data", "raw", "mesn_data.csv"),
                col_types = cols(timestamp = col_character(),
                                 temp = col_double(),
                                 qcflag = col_integer()))

dat <- raw %>%
  filter(!is.na(timestamp), !is.na(temp)) %>%
  mutate(timestamp = ymd_hms(timestamp, tz = "UTC")) %>%
  filter(qcflag == 1) %>%                    # keep only good data
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

range(dat$timestamp)
summary(dat$temp)

# -----------------------------------------------------------------------------
# 2. Daily mean SST
# -----------------------------------------------------------------------------
# NOAA CRW DHW is built from daily SST, so the high-frequency record is
# aggregated to daily means (days with >= 24 samples only).
daily <- dat %>%
  mutate(date = as_date(timestamp)) %>%
  group_by(date) %>%
  summarise(sst = mean(temp, na.rm = TRUE),
            n   = n(), .groups = "drop") %>%
  filter(n >= 24)

ggplot(daily, aes(date, sst)) +
  geom_line(colour = "steelblue") +
  labs(x = NULL, y = "Daily mean SST (°C)",
       title = "MESN buoy daily mean temperature") +
  theme_minimal()

# -----------------------------------------------------------------------------
# 3. Site MMM
# -----------------------------------------------------------------------------
mmm <- 29.9   # site Maximum Monthly Mean (°C)

# -----------------------------------------------------------------------------
# 4. In-situ HotSpot and DHW
# -----------------------------------------------------------------------------
# HotSpot = SST - MMM when SST > MMM (otherwise 0).
# DHW     = sum of HotSpots >= 1 °C over a rolling 84-day window / 7 (°C-weeks).
stress <- daily %>%
  arrange(date) %>%
  mutate(hotspot   = pmax(sst - mmm, 0),
         hotspot_1 = ifelse(hotspot >= 1, hotspot, 0)) %>%
  mutate(dhw = rollapplyr(hotspot_1, width = 84, FUN = sum,
                          fill = NA, partial = TRUE) / 7)

write_csv(stress, here("data", "processed", "mesn_daily_sst_dhw.csv"))

dhw_2024 <- stress %>% filter(year(date) == 2024)

peak_dhw <- dhw_2024 %>% slice_max(dhw, n = 1, with_ties = FALSE)
cat(sprintf("Peak 2024 DHW (in-situ): %.2f °C-weeks on %s\n",
            peak_dhw$dhw, peak_dhw$date))

ggplot(dhw_2024, aes(date, dhw)) +
  geom_area(fill = "firebrick", alpha = 0.4) +
  geom_line(colour = "firebrick") +
  geom_hline(yintercept = c(4, 8), linetype = 2, colour = "grey40") +
  annotate("text", x = min(dhw_2024$date), y = 4.2,
           label = "Bleaching Alert 1 (≥4)", hjust = 0, size = 3) +
  annotate("text", x = min(dhw_2024$date), y = 8.2,
           label = "Bleaching Alert 2 (≥8)", hjust = 0, size = 3) +
  labs(x = NULL, y = "DHW (°C-weeks)",
       title = "2024 accumulated thermal stress (in-situ DHW)") +
  theme_minimal()

# Daily SST against MMM and bleaching threshold (MMM + 1)
ggplot(stress, aes(date, sst)) +
  geom_line(colour = "grey60") +
  geom_hline(yintercept = mmm,     colour = "blue", linetype = 2) +
  geom_hline(yintercept = mmm + 1, colour = "red",  linetype = 2) +
  annotate("text", x = min(stress$date), y = mmm + 0.05,
           label = "MMM (29.9)", hjust = 0, size = 3, colour = "blue") +
  annotate("text", x = min(stress$date), y = mmm + 1.05,
           label = "Bleaching threshold (MMM + 1)", hjust = 0,
           size = 3, colour = "red") +
  labs(x = NULL, y = "Daily mean SST (°C)",
       title = "Daily SST vs. MMM and bleaching threshold") +
  theme_minimal()

# -----------------------------------------------------------------------------
# 5. NOAA CoralTemp 5 km — SST and DHW (2024)
# -----------------------------------------------------------------------------
# CoralTemp v3.1 daily SST (noaacrwsstDaily, field analysed_sst) and DHW
# (noaacrwdhwDaily, field degree_heating_week) for the single 5 km pixel
# nearest the MESN buoy (1.220 °N, 103.843 °E), full year 2024.
#
# The download is slow and depends on the NOAA server, so it runs only if
# data/raw/noaa_coraltemp_2024.csv does not exist; delete that file to
# re-download. Requests are kept small (one pixel, one-month windows, with
# retries) to avoid the ERDDAP "Proxy Error" timeout.
noaa_file <- here("data", "raw", "noaa_coraltemp_2024.csv")

noaa <- read_csv(noaa_file, show_col_types = FALSE) %>%
  select(date, sst, dhw) %>%
  arrange(date)

range(noaa$date)
summary(noaa$sst)

# -----------------------------------------------------------------------------
# 5.1 SST overlay — MESN buoy vs NOAA CoralTemp
# -----------------------------------------------------------------------------
overlay <- bind_rows(
  daily %>% transmute(date, sst, source = "MESN buoy (daily mean)"),
  noaa  %>% transmute(date, sst, source = "NOAA CoralTemp 5 km")
)

ggplot(overlay, aes(date, sst, colour = source)) +
  geom_line(alpha = 0.85, linewidth = 0.4) +
  geom_hline(yintercept = mmm, linetype = 2, colour = "blue") +
  scale_colour_manual(values = c(
    "MESN buoy (daily mean)" = "steelblue",
    "NOAA CoralTemp 5 km"    = "firebrick")) +
  labs(x = NULL, y = "SST (°C)", colour = NULL,
       title = "In-situ vs satellite SST",
       subtitle = sprintf("MMM = %.1f °C (blue)", mmm)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Buoy - NOAA SST agreement on overlapping days
sst_join <- inner_join(
  daily %>% select(date, sst_buoy = sst),
  noaa  %>% select(date, sst_noaa = sst),
  by = "date"
) %>% mutate(diff = sst_buoy - sst_noaa)

cat(sprintf("Buoy - NOAA SST: mean bias %+.2f °C, RMSE %.2f °C, n = %d days\n",
            mean(sst_join$diff, na.rm = TRUE),
            sqrt(mean(sst_join$diff^2, na.rm = TRUE)),
            sum(!is.na(sst_join$diff))))

# -----------------------------------------------------------------------------
# 5.2 Figure: 2024 SST overlay with NOAA DHW inset
# -----------------------------------------------------------------------------
sst_2024 <- bind_rows(
  daily %>% filter(year(date) == 2024) %>%
    transmute(date, sst, source = "MESN buoy (daily mean)"),
  noaa  %>% filter(year(date) == 2024) %>%
    transmute(date, sst, source = "NOAA CoralTemp 5 km")
)

# Main panel: in-situ vs satellite SST
main_sst <- ggplot(sst_2024, aes(date, sst, colour = source)) +
  geom_line(linewidth = 0.5, alpha = 0.9) +
  geom_hline(yintercept = mmm, linetype = 2, colour = "blue") +
  scale_colour_manual(values = c(
    "MESN buoy (daily mean)" = "steelblue",
    "NOAA CoralTemp 5 km"    = "firebrick")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(limits = c(27, 32), breaks = 27:32) +
  labs(x = NULL, y = "SST (°C)", colour = NULL) +
  theme_pubr(base_size = 12) +
  theme(legend.position = "bottom")

# Inset: NOAA CoralTemp DHW (2024)
dhw_inset <- ggplot(noaa %>% filter(year(date) == 2024), aes(date, dhw)) +
  geom_area(fill = "darkorange", alpha = 0.35) +
  geom_line(colour = "darkorange", linewidth = 0.4) +
  geom_hline(yintercept = c(4, 8), linetype = 3, colour = "grey30") +
  annotate("text", x = min(noaa$date), y = 4, vjust = -0.4, hjust = 0,
           label = "Alert 1 (≥4)", size = 2.3, colour = "grey30") +
  annotate("text", x = min(noaa$date), y = 8, vjust = -0.4, hjust = 0,
           label = "Alert 2 (≥8)", size = 2.3, colour = "grey30") +
  scale_x_date(date_labels = "%b") +
  labs(x = NULL, y = "DHW (°C-wk)") +
  theme_pubr(base_size = 7) +
  theme(plot.background = element_rect(fill = "white", colour = "grey70"),
        plot.margin     = margin(3, 4, 3, 3))

# Compose: DHW inset in the empty band below the SST lines (lower-right)
p1 <- main_sst +
  inset_element(dhw_inset,
                left = 0.32, bottom = 0.03, right = 0.85, top = 0.44,
                align_to = "panel")
p1

ggsave(here("output", "figures", "fig_sst_dhw_2024.pdf"), p1,
       width = 9, height = 6)
ggsave(here("output", "figures", "fig_sst_dhw_2024.png"), p1,
       width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 5.3 NOAA CoralTemp DHW (2024)
# -----------------------------------------------------------------------------
peak_noaa <- noaa %>%
  filter(year(date) == 2024) %>%
  slice_max(dhw, n = 1, with_ties = FALSE)

cat(sprintf("Peak 2024 DHW (NOAA CoralTemp): %.2f °C-weeks on %s\n",
            peak_noaa$dhw, peak_noaa$date))

ggplot(noaa %>% filter(year(date) == 2024), aes(date, dhw)) +
  geom_area(fill = "darkorange", alpha = 0.4) +
  geom_line(colour = "darkorange") +
  geom_hline(yintercept = c(4, 8), linetype = 2, colour = "grey40") +
  annotate("text", x = as.Date("2024-01-05"), y = 4.2,
           label = "Bleaching Alert 1 (≥4)", hjust = 0, size = 3) +
  annotate("text", x = as.Date("2024-01-05"), y = 8.2,
           label = "Bleaching Alert 2 (≥8)", hjust = 0, size = 3) +
  labs(x = NULL, y = "DHW (°C-weeks)",
       title = "2024 NOAA CoralTemp DHW (5 km)") +
  theme_minimal()
