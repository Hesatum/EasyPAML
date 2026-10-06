# Tempos do codeml e tempo limite automático

Medição usada para calibrar o tempo limite de cada execução do codeml
(`src/backend/timeouts.py`). O mesmo cálculo vale para a janela e para o CLI.

## Como foi medido

- **codeml 4.9j** (o binário do repositório), um modelo por execução, como o
  EasyPAML roda. Configuração padrão: CodonFreq = 2 (F3×4), ncatG = 10,
  cleandata = 1, ω inicial 0,5, sem warm-start.
- **Dados simulados** com o `evolverNSsites` do PAML 4.10.10: árvores aleatórias de
  10, 30 e 60 táxons; 150, 500 e 1500 códons; modelo M3 com três classes
  (proporções 0,6 / 0,3 / 0,1; ω = 0,05 / 1 / 4), κ = 2. Branch e Branch-site com
  `#1` no primeiro par de folhas irmãs.
- **Máquina**: Intel Core i7-13620H (16 threads), 16 GB, Linux; 8 execuções ao
  mesmo tempo, com `nice 5`. Uma execução sozinha tende a ser um pouco mais rápida.
- 81 execuções (9 tamanhos × 9 modelos); 79 terminaram. As duas de 60 × 1500 com
  M8 e M8a foram cortadas no teto do teste (30 000 s ≈ 8,3 h).

### Tempo medido do codeml (minutos)

| Táxons × códons | M0 | M1a | M2a | M7 | M8 | M8a | Branch | Branch-site | BS nulo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 × 150 | 0.1 | 0.2 | 0.2 | 1.6 | 0.7 | 1.1 | 0.1 | 0.3 | 0.2 |
| 10 × 500 | 0.1 | 0.2 | 0.3 | 1.4 | 1.4 | 2.6 | 0.1 | 0.8 | 1.0 |
| 10 × 1500 | 0.2 | 0.6 | 1.0 | 3.9 | 3.2 | 4.8 | 0.3 | 1.9 | 1.2 |
| 30 × 150 | 0.8 | 2.5 | 3.3 | 11 | 14 | 18 | 1.3 | 5.3 | 5.5 |
| 30 × 500 | 2.1 | 6.5 | 8.8 | 24 | 33 | 33 | 2.2 | 13 | 10 |
| 30 × 1500 | 5.4 | 12 | 22 | 67 | 74 | 82 | 6.6 | 35 | 40 |
| 60 × 150 | 7.8 | 17 | 27 | 57 | 99 | 73 | 6.7 | 61 | 36 |
| 60 × 500 | 20 | 33 | 67 | 203 | 237 | 266 | 16 | 74 | 88 |
| 60 × 1500 | 43 | 84 | 189 | 469 | > 500¹ | > 500¹ | 40 | 240 | 230 |

¹ cortadas no teto do teste (30 000 s); o tempo real é maior.

O tempo cresce muito com o número de táxons e pouco com o comprimento: dobrar os
táxons multiplica o tempo por ~6,6; triplicar os códons, por ~2,2. M7, M8 e M8a
(distribuição beta, 10 categorias) levam de 3,5 a 7 vezes o tempo de M1a/M2a.

## Fórmula

Ajuste por mínimos quadrados em escala log (expoentes comuns a todos os modelos,
uma constante por modelo):

    tempo ≈ T_ref[modelo] × (táxons / 30)^2,73 × (códons / 500)^0,71

| Modelo | M0 | M1a | M2a | M7 | M8 | M8a | Branch | Branch-site | BS nulo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| T_ref (s, 30 táxons × 500 códons) | 140 | 340 | 530 | 1890 | 1950 | 2400 | 160 | 890 | 810 |

**Tempo limite = 5 × a estimativa, no mínimo 30 min.** Modelos fora da tabela usam
o T_ref do Branch-site (se o nome começa com "Branch") ou o maior da tabela.

Folga observada:
- no benchmark, o tempo medido chegou a no máximo 2,4 × a estimativa (M7, 10 × 150);
- em 25 genes reais de *Cereus* (125 execuções de M1a/M2a/M7/M8/M8a, 11 a 24
  táxons, 268 a 1997 códons analisados), o máximo foi 2,4 × a estimativa, ou
  0,47 × o limite; a mediana foi 0,38 × a estimativa.

### Tempo limite automático resultante

| Táxons × códons | M0 | M1a | M2a | M7 | M8 | M8a | Branch | Branch-site | BS nulo |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 × 150 | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min |
| 10 × 500 | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min |
| 10 × 1500 | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min | 30 min |
| 30 × 150 | 30 min | 30 min | 30 min | 1.1 h | 1.2 h | 1.4 h | 30 min | 32 min | 30 min |
| 30 × 500 | 30 min | 30 min | 44 min | 2.6 h | 2.7 h | 3.3 h | 30 min | 1.2 h | 1.1 h |
| 30 × 1500 | 30 min | 1.0 h | 1.6 h | 5.7 h | 5.9 h | 7.3 h | 30 min | 2.7 h | 2.5 h |
| 60 × 150 | 33 min | 1.3 h | 2.1 h | 7.4 h | 7.6 h | 9.4 h | 38 min | 3.5 h | 3.2 h |
| 60 × 500 | 1.3 h | 3.1 h | 4.9 h | 17.4 h | 18.0 h | 22.1 h | 1.5 h | 8.2 h | 7.5 h |
| 60 × 1500 | 2.8 h | 6.8 h | 10.7 h | 38.0 h | 39.2 h | 48.2 h | 3.2 h | 17.9 h | 16.3 h |

## O que protege contra um codeml travado

O tempo limite só corta otimizações que não terminam. Um codeml parado (por
exemplo, esperando um Enter) é encerrado pela **detecção de inatividade**: 5 min
sem usar CPU (`--idle-timeout`, padrão 300 s; 0 desliga).

## Como mudar

- Janela: **Configurações → Tempo limite por modelo (min)**. Vazio = automático.
- CLI: `--timeout SEGUNDOS` (0 = automático, o padrão); em `--config`, a chave
  `"timeout"`.

O limite usado em cada gene × modelo fica no `batch_analysis_log.txt`
(`time limit N s (… taxa × … codons, automatic)`).

## Refazer a medição

Os scripts não ficam no repositório. Em resumo: gerar os alinhamentos com o
`evolverNSsites` (opção 6, `MCcodonNSsites.dat` com os parâmetros acima) e rodar
cada modelo com `CodemlBatchAnalysis` (`n_workers = 1`, `timeout` alto), anotando
o `execution_time` de cada resultado.
