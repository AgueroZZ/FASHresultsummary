


# try fit one model:
index <- 1
datasets[[index]]$offset <- log_size_vec

mod1 <- BayesGP::model_fit(formula = value ~ offset(offset) +
                             f(x, model = "iwp", order = 2,
                               sd.prior = list(h = 1, param = 1)),
                           data = datasets[[index]],
                           family = "poisson")

plot(mod1, select = "x")

mod2 <- BayesGP::model_fit(formula = value ~ offset(offset) +
                             f(x, model = "iwp", order = 2,
                               sd.prior = list(h = 1, param = 1)) +
                             f(Day, model = "iid", sd.prior = list(param = 2)),
                           data = datasets[[index]],
                           family = "poisson")

plot(mod2, select = "x")
