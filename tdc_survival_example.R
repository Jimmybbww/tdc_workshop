library(tidyverse)
library(survival)
library(ggsurvfit)
library(patchwork)

# Functions ====
km_plot = function(x){
  f = reformulate(x, response = 'Surv(time, event)')
  survfit2(f, data = df1) |> 
    ggsurvfit(type = 'risk') +
    #add_risktable(risktable_stats = "{n.risk} ({cum.event})") + 
    add_pvalue() +
    add_confidence_interval() +
    labs(title = x, x = 'Weeks') +
    coord_cartesian(xlim = c(0,60)) +
    scale_y_continuous(limits = c(0, 0.028), breaks = seq(0, 0.03, 0.01)) +
    theme_classic() +
    guides(color  = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.15, 0.8))
  
}
tbl_coxph = function(fit){
  broom::tidy(fit, exponentiate = T, conf.int = T) |> 
    mutate(across(c(estimate,conf.low,conf.high), ~sprintf("%.2f",round(.x,2)))) |> 
    mutate(`HR (95% CI)` = glue::glue("{estimate} ({conf.low}, {conf.high})")) |> 
    select(term, `HR (95% CI)`, p.value)
}

# Read SAS data ====
pth = "https://raw.githubusercontent.com/Jimmybbww/tdc_workshop/master/example.sas7bdat"
df = haven::read_sas(pth)
df = map(df, function(x) {attributes(x)= NULL; x}) |> as.data.frame() |> data.table::setDT()

# Clean data ====
df1 = df |> 
  select(
    id, 
    b_sex_06m, # sex
    starts_with("follow_"), # time
    E6_1_36m, E7_1A_5yf, # event
    e5_2_06m, e3_2_18m, E3_2_36m, E4_2_5yf, # tdc
    b_yy_06m, b_mm_06m, b_dd_06m, # baby age
    m_yy_06m, m_mm_06m, m_dd_06m, # mother age
    f_yy_06m, f_mm_o6m, f_dd_06m, # father age
    ) |> 
  rename(
    sex = b_sex_06m, 
    y_36m = E6_1_36m, 
    gi_6m = e5_2_06m, 
    gi_18m = e3_2_18m, 
    gi_36m = E3_2_36m, 
    gi_60m = E4_2_5yf
    ) |> 
  mutate(
    sex = factor(sex, labels = c('M','F')),
    y_36m = if_else(y_36m=='', '9', y_36m),
    E7_1A_5yf = if_else(E7_1A_5yf=='', '9', E7_1A_5yf),
    y_60m = case_match(E7_1A_5yf, '2'~'1', '1'~'2', .default = E7_1A_5yf),
    bir_b = make_date(as.numeric(b_yy_06m)+1911, b_mm_06m, b_dd_06m),
    bir_m = make_date(as.numeric(m_yy_06m)+1911, m_mm_06m, m_dd_06m),
    bir_f = make_date(as.numeric(f_yy_06m)+1911, f_mm_o6m, f_dd_06m),
    age_mo = bir_m %--% bir_b/years(1),
    age_fa = bir_f %--% bir_b/years(1),
    across(starts_with("gi_"), ~na_if(.x, '9') |> na_if('') ),
    gi = pmax(
      gi_6m, gi_18m, gi_36m, gi_60m, na.rm = T) |> 
      factor(labels = c('No','Yes')),
    across(starts_with("gi_"), ~factor(.x, labels = c('No','Yes')) ),
    follow_max = pmap_dbl(pick(starts_with('follow_')), max, na.rm = T),
    y = paste(y_36m, y_60m, sep = ''),
    time = case_when(
      y %in% c('02','92','20','21','22','29')~ 30/2,
      y %in% c('10','11','12','19','09')~ 30,
      y %in% c('00','01','90','91')~ 54,
      .default = follow_max ), # 99
    event = if_else(y_36m %in% c('1','2') | y_60m %in% c('1','2'), 1, 0)
    ) |> 
  select(
    id, 
    sex, # sex
    starts_with("gi"),
    age_mo, age_fa, # age
    y, time, event
  ) |> 
  filter(time != 0) # time= 0, n= 474

# EDA ====  
skimr::skim(df1)
with(df1, table(y, time))
# table1
table1::table1(~.|sex, data= df1 |> select(-id,-y), render.continuous= "Mean (SD)")

# KM curve ====
p1 = km_plot("gi")
p2 = km_plot("gi_6m")
p3 = km_plot("gi_18m")
p4 = km_plot("gi_36m")
wrap_plots(p1,p2,p3,p4, ncol = 2)

# Cox w/o time-dependent covariates ====
fit0 = coxph(Surv(time, event) ~ gi + sex + age_mo + age_fa, data = df1) 
tbl_coxph(fit0)

# Two-steps for build the time-dependent covariate data ====
## Step 1 ----
base = df1 |> select(id, sex, age_mo, age_fa, time, event, gi_6m)
step1 = tmerge(base, base, id=id, tstop = time)

## Step 2 ----
# Method 1: Using wide-form data directly to build counting process data
step2 = tmerge(
  step1, df1 |> select(id, starts_with('gi')), id = id, 
  gi = tdc(rep(12,nrow(df1)), gi_18m), 
  gi = tdc(rep(30,nrow(df1)), gi_36m)
)

# Method 2: Long-form transformation to build counting process data
tdc = df1 |> 
  select(id, gi_18m, gi_36m) |> 
  pivot_longer(cols = -id, names_to = "time1", values_to = "gi", names_pattern = "(\\d+)") |> 
  mutate(time1 = as.numeric(time1)-6)

step2 = tmerge(step1, tdc, id=id, gi = tdc(time1, gi))

# Cox w/ time-dependent covariates ====
fit = coxph(
  Surv(time, event) ~ gi + gi_6m + sex + age_mo + age_fa, 
  data = step2, cluster = id, ties = 'efron') 
tbl_coxph(fit)



# Highlight functions ====
# `|>` : pipe operator
round(sqrt(sum(1:10)), 2)
1:10 |> sum() |> sqrt() |> round(2)
# dplyr
  # select(), starts_with(), ends_with()
    df1[,1:5]
    df1 |> select(1:5)
    df1 |> select(id, sex, starts_with('age'), ends_with('m'))
    
  # filter()
    df1 |> filter(sex == 'M')
    
  # mutate()
    df1[,c(1,8)] |> mutate(rd.age_mo = round(age_mo))
    
  # rename()
    df1[,c(1,8)] |> rename(age = age_mo)
    
  # if_else()
    df1[,1:2] |> mutate(sex2 = if_else(sex == 'M', 'Male', sex))
    
  # case_when()
    df1[,1:2] |> mutate(sex2 = case_when(sex == 'M'~ 'Male', sex == 'F'~ 'Female', .default = 'Unknown'))

# tidyr
  # pivot_longer()
    df1[,c(1,3:6)] |> pivot_longer(cols = -id, names_to = "time", values_to = "gi")
    # The same as SAS 'proc transpose'
    PROC TRANSPOSE DATA= df1(KEEP = id gi_18m gi_36m) OUT= df1_long(RENAME = (COL1= gi _NAME_= time));
      BY id;
      VAR gi_18m gi_36m; 
    RUN;
    
  # lubridate
    # `%--%` : interval operator
    # make_date()
    x = make_date(2021, 1, 1) %--% today()
    x/years(1); x/months(1); x/weeks(1); x/days(1); x/hours(1); x/minutes(1); x/seconds(1)
  
# survival
  # tmerge()
    
  # coxph()
    m = coxph(Surv(time, status)~ sex, data = lung)
    summary(m)
    
  # survfit
    survfit(m) |> plot()

# ggsurvfit
    km = survfit(Surv(time, status)~ sex, data = lung)
    km
    summary(km)
    plt = ggsurvfit(km) # plotly::ggplotly(plt)
    plt + add_censor_mark() + add_confidence_interval() + add_risktable(risktable_stats = "{n.risk} ({cum.event})")
    
    
    
    
