# Serology for COVID-19
#
# Copyright (c) 2021, Liangliang Wang, Joosung Min, Renny Doig,
#   Lloyd T. Elliott and Caroline Colijn
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Analysis of BC case counts + serology data
# - change in sampling on April 14 -> two separate p_cc parameters

# Preliminaries ----------------------------------------------------------------


rm(list=ls())
library(lubridate)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(rstan)
set.seed(50)


# Data -------------------------------------------------------------------------

# serological data
end_date <- dmy("27-05-2020")
m <- 4
n <- 885

# case counts
raw_data <- read_csv("BCCDC_COVID19_Dashboard_Case_Details.csv")

case_counts <- raw_data %>%
  filter(HA %in% c("Vancouver Coastal", "Fraser")) %>%
  select(Reported_Date) %>%
  group_by(Reported_Date) %>%
  summarise(n()) %>%
  set_names(nm = c("date", "cases")) %>%
  filter(date <= end_date) %>%
  mutate(time = as.numeric(date - min(date)))

start_date <- min(case_counts$date)

case_counts %>%
  ggplot(aes(x=date, y=cases)) + geom_point() +
  geom_vline(xintercept = ymd(start_date + 60))


# Prior settings ---------------------------------------------------------------

# April 14 is t=79 for this dataset

# positive-valued parameters
prior_mu <- c(0.1, 0.1)
prior_eta <- c(0.05, 0.1)
prior_I0 <- c(8, 1)
prior_tau0 <- c(7, 1)
prior_tau <- c(50, 1)

# 0-1-valued parameters
prior_p_cc1 <- c(35, 65) # 0.35 
prior_p_cc2 <- c(65, 35) # 0.65
sens_m <- 0.9
spec_m <- 0.99
spec_v <- 0.009
prior_sens <- c(sens_m*(sens_m*(1-sens_m)/spec_v-1),
                (1-sens_m)*(sens_m*(1-sens_m)/spec_v-1))
prior_spec <- c(spec_m*(spec_m*(1-spec_m)/spec_v-1),
                (1-spec_m)*(spec_m*(1-spec_m)/spec_v-1))


# Fit model --------------------------------------------------------------------

fit_data <- list(T = max(case_counts$time),
                 J = dim(case_counts)[1],
                 time = case_counts$time,
                 m_p = m,
                 n = n,
                 N = 2.85e6,
                 C_cc = case_counts$cases,
                 prior_mu = prior_mu,
                 prior_eta = prior_eta,
                 prior_I0 = prior_I0,
                 prior_tau = prior_tau,
                 prior_tau0 = prior_tau0,
                 prior_p_cc1 = prior_p_cc1,
                 prior_p_cc2 = prior_p_cc2,
                 prior_specificity = prior_spec,
                 prior_sensitivity = prior_sens)

# stan parameters
n_chains = 2
init <- map(1:n_chains, ~ list(I0 = 8,
                               mu = 0.05,
                               tau = 50,
                               eta = 0.01,
                               p_cc1 = 0.35,
                               p_cc2 = 0.65,
                               tau_0 = 7,
                               sensitivity = 0.9,
                               specificity = 0.99))

rstan::rstan_options(auto_write = TRUE)

if(file.exists("fit.rds")){ fit <- readRDS("fit.rds")} else{
  options(mc.cores = 2)
  fit <- rstan::stan(file = "serology_nb.stan",
                      chains = n_chains,
                         init = init,
                         data = fit_data,
                         algorithm = "NUTS",
                         iter = 20000,
                         verbose = T,
                         diagnostic_file = "temp")
  saveRDS(fit, "fit.rds")
}


# Results ----------------------------------------------------------------------

# convergence diagnostic
traceplot(fit, pars=c("mu","eta","I0","tau_0","p_cc1","p_cc2","sensitivity","specificity"))

est_all <- rstan::extract(fit)

est_par <- est_all
est_par[c("I", "lp__")] <- NULL
est_par <- data.frame(est_par) %>%
  gather(key = "Parameter")

est_I <- est_all[["I"]]

est_par %>%
  ggplot(aes(x=value)) +
  geom_histogram(bins=50, colour="grey45") +
  facet_wrap(~Parameter, scale="free") +
  labs(x="", y="") +
  theme(panel.background = element_rect(fill="white", colour="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("BC_histogram.pdf", width=9, height=6, units="in", dpi=900)

cb <- est_par %>%
  group_by(Parameter) %>%
  summarise(mean = round(mean(value),3),
            median = round(quantile(value, probs=0.5),3),
            L95 = round(quantile(value, probs=0.025),3),
            U95 = round(quantile(value, probs=0.975),3))
grid.arrange(tableGrob(cb))
ggsave("credible_intervals.pdf", tableGrob(cb), height=6, width=6, units="in", dpi=900)

# point estimate & credible interval for cumulative counts b/w Jan 27 and May 27
cum_cases <- apply(est_I, 1, sum)
c(mean(cum_cases), quantile(cum_cases, c(0.025, 0.975)))

pt_est <- filter(cb, Parameter%in%c("p_cc1", "p_cc2")) %>%
  bind_rows(data.frame(Parameter="C", mean=mean(cum_cases),
                       median = quantile(cum_cases, 0.5),
                       L95 = quantile(cum_cases, 0.025),
                       U95 = quantile(cum_cases, 0.975)))
rownames(pt_est) <- NULL
grid.arrange(tableGrob(pt_est))
ggsave("point_estimates.pdf", tableGrob(pt_est), height=6, width=6, units="in", dpi=900)

p_ac <- data.frame(est_I) %>%
  set_names(nm=1:fit_data$T) %>%
  gather(key="time", value="I") %>%
  mutate(time = as.numeric(time)) %>%
  group_by(time) %>%
  summarise(mean = mean(I),
            median = quantile(I, probs=0.5),
            L95 = quantile(I, probs=0.025),
            U95 = quantile(I, probs=0.975)) %>%
  mutate(date = ymd(start_date + time - 1)) %>%
  ggplot(aes(x=date, y=mean)) +
  geom_line() +
  geom_ribbon(aes(ymin=L95, ymax=U95), alpha=0.3) +
  xlim(start_date, end_date) +
  labs(x="", y="Infections") +
  theme(panel.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  background_grid(major="xy", colour.major="grey85")

p_cc <- case_counts %>%
  ggplot(aes(x=date, y=cases)) + geom_point() +
  xlim(start_date, end_date) +
  labs(x="Date", y="Observed cases") +
  theme(panel.background = element_rect(fill="white")) +
  background_grid(major="xy", colour.major="grey85")

plot_grid(p_ac, p_cc, align="v", nrow=2, rel_heights=c(2/3,1/3))
ggsave("active_cases.pdf", width=6, height=4, units="in", dpi=900)


priors <- data.frame(prior_mu = prior_mu,
                     prior_eta = prior_eta,
                     prior_I0 = prior_I0,
                     prior_tau = prior_tau,
                     prior_tau0 = prior_tau0,
                     prior_p_cc1 = prior_p_cc1,
                     prior_p_cc2 = prior_p_cc2,
                     prior_specificity = prior_spec,
                     prior_sensitivity = prior_sens)
rownames(priors) <- c("Parameter 1", "Parameter 2")
grid.arrange(tableGrob(priors))

# pairs plot
est_pairs <- rstan::extract(fit, pars=c("I0","p_cc1","p_cc2","mu"))
pdf("pairs", width=9, height=8)
pairs(data.frame(est_pairs), cex=0.1)
dev.off()
