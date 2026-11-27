# AGENTS.md

## Contexto do projeto

Trabalho prático da disciplina **Métodos Estatísticos de Previsão**, escrito em português e desenvolvido em R/Quarto. O objetivo é comparar métodos de previsão para uma série mensal de 120 observações, de janeiro de 2010 a dezembro de 2019.

O relatório apresenta a variável como número mensal de casos de suicídio. O CSV, porém, contém somente data (`index`) e contagem (`n`): a fonte primária, a abrangência geográfica e a definição do caso ainda precisam ser documentadas. Não invente esses dados nem conclusões substantivas que dependam deles.

## Arquivos

- `Trabalho Prático.qmd`: relatório e análise principal; é a fonte do PDF.
- `total_por_ano_e_mes.csv`: dados locais usados na análise.
- `box_jenkins_arima.R`: script separado com a metodologia Box-Jenkins; use-o como referência para a seção ARIMA do relatório.
- `AulaPratica1.R` e `AulaPratica2.R`: exemplos das aulas práticas.
- `Trabalho-Prático.pdf`: PDF presente no diretório; confirme se corresponde à versão atual do `.qmd`.

## Estrutura da análise

- **Dados:** série mensal; treino de janeiro de 2010 a junho de 2018 (102 meses) e teste de julho de 2018 a dezembro de 2019 (18 meses). Preserve a ordem cronológica e use a mesma divisão para todos os modelos.
- **Regressão:** compare tendência linear, tendência com sazonalidade e tendência quadrática com sazonalidade. Examine ajuste e resíduos; escolha o representante no treino com AICc e diagnóstico.
- **Suavização exponencial:** compare suavização simples, Holt e Holt-Winters aditivo. Examine ajuste e resíduos; escolha o representante no treino com AICc e diagnóstico.
- **ARIMA:** siga Box-Jenkins. Faça testes de estacionariedade e decida as diferenças sazonais e regulares no treino; examine ACF e PACF; estime candidatos; diagnostique os resíduos e revise a identificação quando necessário. Entre candidatos adequados, use o AIC para escolher o representante.
- **Comparação final:** avalie os três representantes nos mesmos 18 meses de teste com MAE, RMSE, MASE e gráficos das previsões. Considere também a qualidade dos resíduos.

## Relatório e verificação

- Deixe o código executável no `.qmd`, mas oculto no PDF. Mostre apenas resultados, tabelas e figuras necessários, com interpretação breve.
- Use títulos, eixos e legendas em português. Evite tabelas com nomes de linha duplicados e formate valores-p sem apresentar `0,000` como probabilidade exata.
- Não use o teste para ajustar ou escolher candidatos dentro das famílias. Não altere os dados brutos sem explicar a transformação. Use caminhos relativos ao repositório.
- Antes da entrega, renderize o `.qmd` completo e confira cálculos, avisos, paginação, legendas, tabelas e a posição das figuras. Os números citados no texto devem corresponder à saída gerada.

## Pendências conhecidas

- Documentar a fonte, a definição da variável e a abrangência geográfica.
- Redigir a conclusão, que está vazia no `.qmd` atual.
- Trocar o caminho absoluto do CSV em `box_jenkins_arima.R` por um caminho relativo antes de usar o script fora desta máquina.
