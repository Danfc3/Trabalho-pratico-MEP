# Trabalho Prático de Métodos Estatísticos de Previsão
# Código completo da análise do relatório (Trabalho Prático.qmd).
# Executar na pasta que contém o arquivo total_por_ano_e_mes.csv.

## ---- setup ----
library(dplyr)
library(lubridate)
library(forecast)
library(ggplot2)
library(tseries)
library(knitr)

meses <- c("Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
           "Jul", "Ago", "Set", "Out", "Nov", "Dez")
meses_ext <- c("janeiro", "fevereiro", "março", "abril", "maio", "junho",
               "julho", "agosto", "setembro", "outubro", "novembro",
               "dezembro")

## formatação dos números no padrão brasileiro
f <- function(x, d = 2) formatC(x, format = "f", digits = d, decimal.mark = ",")
## vírgula decimal dentro de fórmulas
fm <- function(x, d = 2) gsub(",", "{,}", f(x, d))
fp <- function(p) ifelse(p < 0.001, "< 0,001", f(p, 3))
pv <- function(p) {
  ifelse(p < 0.001, "p-valor < 0,001", paste("p-valor =", f(p, 3)))
}
fp2 <- function(p1, p2) {
  if (p1 < 0.001 && p2 < 0.001) "abaixo de 0,001 nos dois casos"
  else paste(fp(p1), "e", fp(p2))
}
fp_adf <- function(p) {
  ifelse(p <= 0.01, "< 0,01", ifelse(p >= 0.99, "> 0,99", f(p, 3)))
}
fp_kpss <- function(p) {
  ifelse(p <= 0.01, "< 0,01", ifelse(p >= 0.1, "> 0,10", f(p, 3)))
}
tabela <- function(x, digitos = 2, ...) {
  knitr::kable(x, digits = digitos, row.names = FALSE,
               format.args = list(decimal.mark = ","), ...)
}

## testes nos resíduos: Ljung-Box com 24 defasagens (gl = parâmetros ARMA)
lb <- function(e, gl = 0) {
  Box.test(e, lag = 24, type = "Ljung-Box", fitdf = gl)$p.value
}
sw <- function(e) shapiro.test(as.numeric(e))$p.value
n_fora <- function(e) {
  sum(abs(Acf(e, lag.max = 24, plot = FALSE)$acf[-1]) > 1.96 / sqrt(length(e)))
}
frase_fora <- function(e) {
  k <- n_fora(e)
  if (k == 0) "nenhuma das 24 autocorrelações fica fora dos limites"
  else if (k == 1) "só 1 das 24 autocorrelações fica fora dos limites"
  else if (k == 24) "todas as 24 autocorrelações ficam fora dos limites"
  else paste(k, "das 24 autocorrelações ficam fora dos limites")
}

## resíduos do ARIMA sem as d + 12D primeiras observações, que ficam
## praticamente nulas por causa da inicialização das diferenças
residuos_arima <- function(modelo) {
  o <- arimaorder(modelo)
  k <- o[2] + if (length(o) > 3) 12 * o[5] else 0
  e <- residuals(modelo)
  window(e, start = time(e)[k + 1])
}
nome_arima <- function(modelo) {
  o <- arimaorder(modelo)
  if (length(o) == 3 || all(o[4:6] == 0)) {
    sprintf("ARIMA(%d,%d,%d)", o[1], o[2], o[3])
  } else {
    sprintf("ARIMA(%d,%d,%d)(%d,%d,%d)[12]", o[1], o[2], o[3], o[4], o[5], o[6])
  }
}
D_de <- function(modelo) {
  o <- arimaorder(modelo)
  if (length(o) == 3) 0 else o[5]
}

## medidas na base de teste e cobertura do intervalo de 95%
medidas <- function(prev) {
  accuracy(prev, serie_teste)[2, c("RMSE", "MAE", "MAPE", "MASE")]
}
n_dentro <- function(prev) {
  sum(serie_teste >= prev$lower[, 2] & serie_teste <= prev$upper[, 2])
}
cobertura <- function(prev) 100 * n_dentro(prev) / length(serie_teste)
frase_dentro <- function(prev) {
  k <- n_dentro(prev)
  if (k == length(serie_teste)) {
    paste("todos os", k, "valores do teste ficam")
  } else {
    paste(k, "dos", length(serie_teste), "valores do teste ficam")
  }
}

## força da sazonalidade, mesmo critério do nsdiffs (diferenciar se > 0,64)
forca_sazonal <- function(x) {
  decomp <- mstl(x)
  resto <- decomp[, "Remainder"]
  sazonal <- decomp[, grep("Season", colnames(decomp))]
  max(0, min(1, 1 - var(resto) / var(resto + sazonal)))
}

## painel de diagnóstico dos resíduos, seguindo o roteiro da aula:
## resíduos no tempo, ACF, PACF, normalidade e Ljung-Box
diagnostico <- function(e, gl = 0) {
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 3.5, 1))
  plot(e, xlab = "Ano", ylab = "Resíduo", main = "Resíduos no tempo")
  abline(h = 0, lty = 2, col = "grey40")
  Acf(e, lag.max = 24, main = "ACF dos resíduos")
  Pacf(e, lag.max = 24, main = "PACF dos resíduos")
  hist(as.numeric(e), breaks = 15, freq = FALSE, main = "Histograma",
       xlab = "Resíduo", ylab = "Densidade")
  curve(dnorm(x, mean(e), sd(e)), add = TRUE, col = "red")
  qqnorm(as.numeric(e), main = "Gráfico Q-Q", xlab = "Quantis teóricos",
         ylab = "Quantis dos resíduos")
  qqline(as.numeric(e), col = "red")
  defasagens <- (gl + 1):24
  p_valores <- sapply(defasagens, function(k) {
    Box.test(e, lag = k, type = "Ljung-Box", fitdf = gl)$p.value
  })
  plot(defasagens, p_valores, ylim = c(0, 1), pch = 19, xlab = "Defasagem",
       ylab = "p-valor", main = "Ljung-Box")
  abline(h = 0.05, lty = 2, col = "blue")
  par(op)
}

## ---- dados ----
desviv <- read.csv("total_por_ano_e_mes.csv")
names(desviv)[names(desviv) == "index"] <- "Data"
names(desviv)[names(desviv) == "n"] <- "Casos"
desviv <- desviv %>%
  mutate(
    Data = as.Date(Data),
    Ano = year(Data),
    Mes = month(Data)
  )

ts_des <- ts(desviv$Casos, start = c(2010, 1), frequency = 12)

## verificação: meses consecutivos, sem datas repetidas e sem valores ausentes
datas_esperadas <- seq(as.Date("2010-01-01"), as.Date("2019-12-01"),
                       by = "month")
stopifnot(identical(desviv$Data, datas_esperadas),
          sum(is.na(desviv$Casos)) == 0,
          sum(duplicated(desviv$Data)) == 0)

## valores atípicos pela distância interquartílica
quartis <- quantile(desviv$Casos, c(0.25, 0.75))
aiq <- diff(quartis)
limites_aiq <- c(quartis[1] - 1.5 * aiq, quartis[2] + 1.5 * aiq)
n_atipicos <- sum(desviv$Casos < limites_aiq[1] | desviv$Casos > limites_aiq[2])

## ---- fig-serie ----
autoplot(ts_des) +
  labs(x = "Ano", y = "Número de casos") +
  theme_minimal()

## ---- tbl-descritivas ----
resumo <- data.frame(
  Medida = c("Mínimo", "1º quartil", "Mediana", "Média", "3º quartil",
             "Máximo", "Desvio padrão"),
  Valor = c(min(ts_des), quantile(ts_des, 0.25), median(ts_des), mean(ts_des),
            quantile(ts_des, 0.75), max(ts_des), sd(ts_des))
)
tabela(resumo, digitos = 1)

## ---- tbl-anual ----
anual <- desviv %>%
  group_by(Ano) %>%
  summarise(Total = sum(Casos), Media = mean(Casos), DP = sd(Casos)) %>%
  mutate(Ano = as.character(Ano))
tabela(anual, digitos = 1,
       col.names = c("Ano", "Total", "Média mensal", "Desvio padrão"))

## ---- fig-decomposicao ----
decomposicao <- decompose(ts_des, type = "additive")
resto <- na.omit(decomposicao$random)
z_resto <- resto / sd(resto)
pos_max <- which.max(abs(z_resto))
mes_max <- paste(meses_ext[cycle(z_resto)[pos_max]], "de",
                 floor(time(z_resto)[pos_max] + 0.001))

op <- par(mfrow = c(4, 1), mar = c(2, 4, 2, 1))
plot(decomposicao$x, main = "Série observada", xlab = "", ylab = "Casos")
plot(decomposicao$trend, main = "Tendência", xlab = "", ylab = "Casos")
plot(decomposicao$seasonal, main = "Sazonalidade", xlab = "", ylab = "Casos")
plot(decomposicao$random, main = "Componente irregular", xlab = "Ano",
     ylab = "Casos")
par(op)

## ---- tbl-mensal ----
mensal <- desviv %>%
  group_by(Mes) %>%
  summarise(Media = mean(Casos)) %>%
  mutate(Mes = meses, Efeito = decomposicao$figure)
tabela(mensal, digitos = 1, col.names = c("Mês", "Média", "Efeito sazonal"))

## ---- fig-acf-pacf ----
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3.5, 1))
Acf(ts_des, lag.max = 24, main = "ACF")
Pacf(ts_des, lag.max = 24, main = "PACF")
par(op)
acf_orig <- Acf(ts_des, lag.max = 24, plot = FALSE)$acf[-1]
pacf_orig <- Pacf(ts_des, lag.max = 24, plot = FALSE)$acf

## ---- tbl-testes ----
ts_des_diff <- diff(ts_des, differences = 1)
adf_orig <- tseries::adf.test(ts_des)
kpss_orig <- tseries::kpss.test(ts_des)
adf_dif <- tseries::adf.test(ts_des_diff)
kpss_dif <- tseries::kpss.test(ts_des_diff)

testes <- data.frame(
  Serie = c("Original", "Original", "1ª diferença", "1ª diferença"),
  Teste = rep(c("ADF (H0: raiz unitária)", "KPSS (H0: estacionária)"), 2),
  Estatistica = c(adf_orig$statistic, kpss_orig$statistic,
                  adf_dif$statistic, kpss_dif$statistic),
  p = c(fp_adf(adf_orig$p.value), fp_kpss(kpss_orig$p.value),
        fp_adf(adf_dif$p.value), fp_kpss(kpss_dif$p.value))
)
tabela(testes, digitos = 3,
       col.names = c("Série", "Teste", "Estatística", "p-valor"))

## ---- fig-diferenca ----
ggtsdisplay(ts_des_diff, lag.max = 36, main = "", xlab = "Ano",
            ylab = "Diferença", theme = theme_minimal())

## ---- divisao ----
n <- length(ts_des)
n_teste <- min(18, round(0.2 * n))
serie_treino <- window(ts_des, end = time(ts_des)[n - n_teste])
serie_teste <- window(ts_des, start = time(ts_des)[n - n_teste + 1])

## ---- fig-treino-teste ----
autoplot(serie_treino, series = "Treino") +
  autolayer(serie_teste, series = "Teste") +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- regressao ----
modelo_tendencia <- tslm(serie_treino ~ trend)
modelo_reg_sazonal <- tslm(serie_treino ~ trend + season)
modelo_reg_quadratico <- tslm(serie_treino ~ trend + I(trend^2) + season)

regressoes <- list(
  "Tendência linear" = modelo_tendencia,
  "Tendência + sazonalidade" = modelo_reg_sazonal,
  "Tendência quadrática + sazonalidade" = modelo_reg_quadratico
)
tab_reg <- data.frame(
  Modelo = names(regressoes),
  R2 = sapply(regressoes, function(m) summary(m)$adj.r.squared),
  AICc = sapply(regressoes, function(m) CV(m)["AICc"]),
  LB = sapply(regressoes, function(m) lb(residuals(m))),
  SW = sapply(regressoes, function(m) sw(residuals(m)))
)
acf12_tendencia <- Acf(residuals(modelo_tendencia), plot = FALSE)$acf[13]

## ---- fig-reg-residuos ----
op <- par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))
for (nome in names(regressoes)) {
  e <- residuals(regressoes[[nome]])
  plot(e, main = nome, xlab = "", ylab = "Resíduo")
  abline(h = 0, lty = 2, col = "grey40")
  lines(lowess(time(e), e), col = "red", lwd = 2)
}
par(op)

## ---- tbl-regressao ----
tab_reg_saida <- tab_reg
tab_reg_saida$LB <- fp(tab_reg$LB)
tab_reg_saida$SW <- fp(tab_reg$SW)
tabela(tab_reg_saida, digitos = c(0, 3, 2, 0, 0),
       col.names = c("Modelo", "$R^2$ ajustado", "AICc", "p-valor Ljung-Box",
                     "p-valor Shapiro-Wilk"))

## ---- tbl-reg-coef ----
## restrição de soma zero nos efeitos mensais: mesmo ajuste, outra leitura
sazonalidade <- factor(cycle(serie_treino))
contrasts(sazonalidade) <- contr.sum(12)
modelo_reg_quadratico_sz <- tslm(serie_treino ~ trend + I(trend^2) +
                                   sazonalidade)
stopifnot(isTRUE(all.equal(as.numeric(fitted(modelo_reg_quadratico)),
                           as.numeric(fitted(modelo_reg_quadratico_sz)))))

coefs <- summary(modelo_reg_quadratico_sz)$coefficients
idx <- grep("sazonalidade", rownames(coefs))
## efeito de dezembro = menos a soma dos outros 11 efeitos
dez_est <- -sum(coefs[idx, 1])
dez_ep <- sqrt(sum(vcov(modelo_reg_quadratico_sz)[idx, idx]))
dez_p <- 2 * pt(-abs(dez_est / dez_ep),
                df = modelo_reg_quadratico_sz$df.residual)
efeitos <- c(coefs[idx, 1], dez_est)

tab_coef <- data.frame(
  Termo = c("Intercepto", "$t$", "$t^2$", paste("Efeito de", meses)),
  Estimativa = c(coefs[1:3, 1], efeitos),
  EP = c(coefs[1:3, 2], coefs[idx, 2], dez_ep),
  p = fp(c(coefs[1:3, 4], coefs[idx, 4], dez_p))
)
tabela(tab_coef, digitos = c(0, 3, 3, 0),
       col.names = c("Termo", "Estimativa", "Erro padrão", "p-valor"))

## ---- fig-reg-ajuste ----
autoplot(serie_treino, series = "Observado") +
  autolayer(fitted(modelo_reg_quadratico),
            series = "Tendência quadrática + sazonalidade") +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- fig-reg-diagnostico ----
diagnostico(residuals(modelo_reg_quadratico))

## ---- prev-regressao ----
prev_tendencia <- forecast(modelo_tendencia, h = n_teste)
prev_sazonal <- forecast(modelo_reg_sazonal, h = n_teste)
prev_quadratico <- forecast(modelo_reg_quadratico, h = n_teste)

prev_reg <- list(
  "Tendência linear" = prev_tendencia,
  "Tendência + sazonalidade" = prev_sazonal,
  "Tendência quadrática + sazonalidade" = prev_quadratico
)
tab_reg_teste <- data.frame(Modelo = names(prev_reg),
                            t(sapply(prev_reg, medidas)),
                            Cobertura = sapply(prev_reg, cobertura))
me_linear <- accuracy(prev_tendencia, serie_teste)[2, "ME"]
me_sazonal <- accuracy(prev_sazonal, serie_teste)[2, "ME"]

## ---- fig-reg-previsao ----
autoplot(window(ts_des, start = c(2016, 1)), series = "Observado") +
  geom_vline(xintercept = time(serie_teste)[1], linetype = "dashed",
             colour = "grey40") +
  autolayer(prev_tendencia$mean, series = "Tendência linear") +
  autolayer(prev_sazonal$mean, series = "Tendência + sazonalidade") +
  autolayer(prev_quadratico$mean,
            series = "Tendência quadrática + sazonalidade") +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- tbl-reg-teste ----
tabela(tab_reg_teste, digitos = 2,
       col.names = c("Modelo", "RMSE", "MAE", "MAPE (%)", "MASE",
                     "Cobertura IP 95% (%)"))

## ---- suavizacao ----
modelo_ANN <- ets(serie_treino, model = "ANN", damped = FALSE)
modelo_AAN <- ets(serie_treino, model = "AAN", damped = FALSE)
modelo_AAA <- ets(serie_treino, model = "AAA", damped = FALSE)

suavizacoes <- list(
  "Suavização simples, ETS(A,N,N)" = modelo_ANN,
  "Holt, ETS(A,A,N)" = modelo_AAN,
  "Holt-Winters aditivo, ETS(A,A,A)" = modelo_AAA
)
tab_ets <- data.frame(
  Modelo = names(suavizacoes),
  AICc = sapply(suavizacoes, function(m) m$aicc),
  LB = sapply(suavizacoes, function(m) lb(residuals(m))),
  SW = sapply(suavizacoes, function(m) sw(residuals(m)))
)

## ---- tbl-suavizacao ----
tab_ets_saida <- tab_ets
tab_ets_saida$LB <- fp(tab_ets$LB)
tab_ets_saida$SW <- fp(tab_ets$SW)
tabela(tab_ets_saida, digitos = c(0, 2, 0, 0),
       col.names = c("Modelo", "AICc", "p-valor Ljung-Box",
                     "p-valor Shapiro-Wilk"))

## ---- tbl-ets-par ----
par_AAA <- modelo_AAA$par
tab_par <- data.frame(
  Parametro = c("$\\alpha$ (nível)", "$\\beta$ (tendência)",
                "$\\gamma$ (sazonalidade)"),
  Estimativa = as.numeric(par_AAA[c("alpha", "beta", "gamma")])
)
tabela(tab_par, digitos = 4, col.names = c("Parâmetro", "Estimativa"))

## ---- fig-ets-ajuste ----
autoplot(serie_treino, series = "Observado") +
  autolayer(fitted(modelo_AAA), series = "Holt-Winters") +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- fig-ets-diagnostico ----
diagnostico(residuals(modelo_AAA))

## ---- prev-suavizacao ----
prev_ANN <- forecast(modelo_ANN, h = n_teste)
prev_AAN <- forecast(modelo_AAN, h = n_teste)
prev_AAA <- forecast(modelo_AAA, h = n_teste)

prev_ets <- list(
  "Suavização simples" = prev_ANN,
  "Holt" = prev_AAN,
  "Holt-Winters" = prev_AAA
)
tab_ets_teste <- data.frame(Modelo = names(prev_ets),
                            t(sapply(prev_ets, medidas)),
                            Cobertura = sapply(prev_ets, cobertura))

## ---- fig-ets-previsao ----
autoplot(window(ts_des, start = c(2016, 1)), series = "Observado") +
  geom_vline(xintercept = time(serie_teste)[1], linetype = "dashed",
             colour = "grey40") +
  autolayer(prev_ANN$mean, series = "Suavização simples") +
  autolayer(prev_AAN$mean, series = "Holt") +
  autolayer(prev_AAA$mean, series = "Holt-Winters") +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- tbl-ets-teste ----
tabela(tab_ets_teste, digitos = 2,
       col.names = c("Modelo", "RMSE", "MAE", "MAPE (%)", "MASE",
                     "Cobertura IP 95% (%)"))

## ---- arima-identificacao ----
serie_treino_diff <- diff(serie_treino, differences = 1)
adf_treino <- tseries::adf.test(serie_treino)
kpss_treino <- tseries::kpss.test(serie_treino)
adf_treino_dif <- tseries::adf.test(serie_treino_diff)
kpss_treino_dif <- tseries::kpss.test(serie_treino_diff)
fs_treino <- forca_sazonal(serie_treino)
D_forca <- nsdiffs(serie_treino)
D_ocsb <- nsdiffs(serie_treino, test = "ocsb")
acf_dif <- Acf(serie_treino_diff, lag.max = 36, plot = FALSE)$acf[-1]
pacf_dif <- Pacf(serie_treino_diff, lag.max = 36, plot = FALSE)$acf
lim_dif <- 1.96 / sqrt(length(serie_treino_diff))

## ---- tbl-arima-testes ----
tab_ident <- data.frame(
  Teste = c("ADF na série", "KPSS na série", "ADF na 1ª diferença",
            "KPSS na 1ª diferença",
            "Força da sazonalidade (diferenciar se > 0,64)",
            "OCSB (diferença sazonal)"),
  Resultado = c(paste("p-valor", fp_adf(adf_treino$p.value)),
                paste("p-valor", fp_kpss(kpss_treino$p.value)),
                paste("p-valor", fp_adf(adf_treino_dif$p.value)),
                paste("p-valor", fp_kpss(kpss_treino_dif$p.value)),
                f(fs_treino, 3),
                ifelse(D_ocsb == 0, "não indica diferença",
                       "indica diferença")),
  Indicacao = c("rejeita raiz unitária (teste com tendência)",
                "não estacionária", "estacionária", "estacionária: d = 1",
                paste("D =", D_forca), paste("D =", D_ocsb))
)
tabela(tab_ident, col.names = c("Teste", "Resultado", "Indicação"))

## ---- fig-arima-dif ----
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3.5, 1))
Acf(serie_treino_diff, lag.max = 36, main = "ACF da 1ª diferença")
Pacf(serie_treino_diff, lag.max = 36, main = "PACF da 1ª diferença")
par(op)

## ---- arima-rodada1 ----
modelo1 <- Arima(serie_treino, order = c(0, 1, 1))  ## ARIMA(0,1,1)
modelo2 <- Arima(serie_treino, order = c(1, 1, 0))  ## ARIMA(1,1,0)
modelo_auto <- auto.arima(serie_treino,
                          d = 1,
                          max.p = 3,
                          max.q = 3,
                          seasonal = F,
                          allowdrift = F)
residuos1 <- residuos_arima(modelo1)
acf_res1 <- Acf(residuos1, lag.max = 36, plot = FALSE)$acf[-1]

## ---- fig-arima-res1 ----
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3.5, 1))
Acf(residuos1, lag.max = 36, main = "ACF dos resíduos")
Pacf(residuos1, lag.max = 36, main = "PACF dos resíduos")
par(op)

## ---- arima-residuos ----
res_sar <- Arima(residuos1, seasonal = c(1, 0, 0), include.mean = FALSE)
res_sma <- Arima(residuos1, seasonal = c(0, 1, 1))
## parte MA com o sinal da disciplina: a_t - Theta a_(t-12)
tab_res <- data.frame(
  Modelo = c("ARIMA(0,0,0)(1,0,0)[12]", "ARIMA(0,0,0)(0,1,1)[12]"),
  Parametro = c("$\\Phi_1$ (AR sazonal)", "$\\Theta_1$ (MA sazonal)"),
  Estimativa = c(coef(res_sar), -coef(res_sma)),
  LB = fp(c(lb(residuos_arima(res_sar), 1), lb(residuos_arima(res_sma), 1)))
)

## ---- tbl-arima-residuos ----
tabela(tab_res, digitos = c(0, 0, 3, 0),
       col.names = c("Modelo para os resíduos", "Parâmetro", "Estimativa",
                     "p-valor Ljung-Box"))

## ---- fig-arima-dif2 ----
serie_treino_dif2 <- diff(serie_treino_diff, lag = 12)
acf_dif2 <- Acf(serie_treino_dif2, lag.max = 36, plot = FALSE)$acf[-1]
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3.5, 1))
Acf(serie_treino_dif2, lag.max = 36, main = "ACF com as duas diferenças")
Pacf(serie_treino_dif2, lag.max = 36, main = "PACF com as duas diferenças")
par(op)

## ---- arima-candidatos ----
modelo3 <- Arima(serie_treino, order = c(0, 1, 1), seasonal = c(1, 0, 0))
modelo4 <- Arima(serie_treino, order = c(1, 1, 0), seasonal = c(1, 0, 0))
modelo5 <- Arima(serie_treino, order = c(0, 1, 1), seasonal = c(0, 1, 1))
modelo6 <- Arima(serie_treino, order = c(1, 1, 0), seasonal = c(0, 1, 1))
modelo_auto_saz <- auto.arima(serie_treino,
                              d = 1,
                              max.p = 3,
                              max.q = 3,
                              seasonal = T,
                              allowdrift = F)

arimas <- list(modelo1, modelo2, modelo3, modelo4, modelo5, modelo6)
names(arimas) <- sapply(arimas, nome_arima)
gl_arima <- sapply(arimas, function(m) length(coef(m)))
tab_arima <- data.frame(
  Modelo = names(arimas),
  D = sapply(arimas, D_de),
  AICc = sapply(arimas, function(m) m$aicc),
  LB = mapply(function(m, k) lb(residuos_arima(m), k), arimas, gl_arima)
)
tab_arima$Adequado <- ifelse(tab_arima$LB >= 0.05, "sim", "não")

## ---- tbl-arima ----
tab_arima_saida <- tab_arima
tab_arima_saida$LB <- fp(tab_arima$LB)
tabela(tab_arima_saida, digitos = c(0, 0, 2, 0, 0),
       col.names = c("Modelo", "D", "AICc", "p-valor Ljung-Box", "Adequado"))

## ---- tbl-arima-coef ----
modelo_sarima <- modelo5
## o R escreve a_t + theta a_(t-1); a disciplina, a_t - theta a_(t-1)
est <- -coef(modelo_sarima)
ep <- sqrt(diag(modelo_sarima$var.coef))
tab_sarima <- data.frame(
  Parametro = c("$\\theta_1$ (MA)", "$\\Theta_1$ (MA sazonal)"),
  Estimativa = est,
  EP = ep,
  Li = est - 1.96 * ep,
  Ls = est + 1.96 * ep,
  p = fp(2 * pnorm(-abs(est / ep)))
)
tabela(tab_sarima, digitos = c(0, 3, 3, 3, 3, 0),
       col.names = c("Parâmetro", "Estimativa", "Erro padrão",
                     "IC 95%: limite inferior", "IC 95%: limite superior",
                     "p-valor"))
residuos_sarima <- residuos_arima(modelo_sarima)

## ---- fig-arima-diagnostico ----
diagnostico(residuos_sarima, gl = 2)

## ---- prev-arima ----
prev_arima <- lapply(arimas, forecast, h = n_teste)
prev_sarima <- prev_arima[["ARIMA(0,1,1)(0,1,1)[12]"]]
tab_arima_teste <- data.frame(Modelo = names(prev_arima),
                              t(sapply(prev_arima, medidas)),
                              Cobertura = sapply(prev_arima, cobertura))

## ---- fig-arima-previsao ----
autoplot(window(ts_des, start = c(2016, 1)), series = "Observado") +
  geom_vline(xintercept = time(serie_teste)[1], linetype = "dashed",
             colour = "grey40") +
  autolayer(prev_arima[["ARIMA(0,1,1)"]]$mean, series = "ARIMA(0,1,1)") +
  autolayer(prev_arima[["ARIMA(1,1,0)(1,0,0)[12]"]]$mean,
            series = "ARIMA(1,1,0)(1,0,0)[12]") +
  autolayer(prev_sarima$mean, series = "ARIMA(0,1,1)(0,1,1)[12]") +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- tbl-arima-teste ----
tabela(tab_arima_teste, digitos = 2,
       col.names = c("Modelo", "RMSE", "MAE", "MAPE (%)", "MASE",
                     "Cobertura IP 95% (%)"))

## ---- comparacao ----
prev_ingenuo <- snaive(serie_treino, h = n_teste)
representantes <- list(
  "Regressão: tendência quadrática + sazonalidade" = prev_quadratico,
  "Suavização: Holt-Winters aditivo" = prev_AAA,
  "ARIMA(0,1,1)(0,1,1)[12]" = prev_sarima,
  "Referência: ingênuo sazonal" = prev_ingenuo
)
## RMSE no treino (para o ARIMA, sem os resíduos da inicialização)
rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
rmse_treino <- c(rmse(residuals(modelo_reg_quadratico)),
                 rmse(residuals(modelo_AAA)),
                 rmse(residuos_sarima),
                 rmse(residuals(prev_ingenuo)))
comparacao <- data.frame(Modelo = names(representantes),
                         RMSE_treino = rmse_treino,
                         t(sapply(representantes, medidas)),
                         Cobertura = sapply(representantes, cobertura))

## ---- tbl-final ----
tabela(comparacao, digitos = 2,
       col.names = c("Modelo", "RMSE no treino", "RMSE", "MAE", "MAPE (%)",
                     "MASE", "Cobertura IP 95% (%)"))

## ---- fig-final ----
autoplot(window(ts_des, start = c(2016, 1)), series = "Observado") +
  geom_vline(xintercept = time(serie_teste)[1], linetype = "dashed",
             colour = "grey40") +
  autolayer(prev_quadratico$mean,
            series = "Regressão quadrática + sazonalidade") +
  autolayer(prev_AAA$mean, series = "Holt-Winters") +
  autolayer(prev_sarima$mean, series = "ARIMA(0,1,1)(0,1,1)[12]") +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

## ---- tbl-psi ----
th <- coef(modelo_sarima)  # notação do R
psi <- ARMAtoMA(ar = c(1, rep(0, 10), 1, -1),
                ma = c(th[1], rep(0, 10), th[2], th[1] * th[2]),
                lag.max = n_teste - 1)
V_h <- modelo_sarima$sigma2 * cumsum(c(1, psi^2))
prev_sarima95 <- forecast(modelo_sarima, h = n_teste, level = 95)
hs <- c(1, 2, 3, 12, 13, 18)
tab_psi <- data.frame(
  h = hs,
  psi = c(1, psi)[hs],
  V = V_h[hs],
  Formula = 1.96 * sqrt(V_h[hs]),
  R = as.numeric(prev_sarima95$upper - prev_sarima95$mean)[hs]
)
tabela(tab_psi, digitos = c(0, 3, 1, 2, 2),
       col.names = c("h", "$\\psi_{h-1}$", "$V(h)$",
                     "Semiamplitude pela fórmula", "Semiamplitude no R"))
dif_psi <- max(abs(tab_psi$R / tab_psi$Formula - 1))

## ---- fig-intervalo ----
autoplot(prev_quadratico, include = 30) +
  autolayer(serie_teste, series = "Observado no teste") +
  labs(title = NULL, x = "Ano", y = "Número de casos", colour = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

