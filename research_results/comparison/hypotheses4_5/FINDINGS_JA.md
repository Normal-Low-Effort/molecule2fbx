# 仮説4・5 既存ensemble解析

新規ORCA計算は行っていない。主要集計は共通骨格RMSD clusterごとの最低電子エネルギー代表を1個だけ採用し、298.15 Kで電子エネルギー重み付けした感度解析である。これは厳密な自由エネルギーpopulationではない。

## 仮説4：配座集団効果

- Bz: 共通骨格 9 cluster、反応中心 4 cluster、5 kJ/mol以内 4構造、電子Eによる実効配座数 4.55。
- SB: 共通骨格 9 cluster、反応中心 4 cluster、5 kJ/mol以内 4構造、電子Eによる実効配座数 4.58。
- SBの全重原子cluster増加を、そのままLSD–benzoyl骨格の多様性増加とは扱えない。TMS末端を除くと両者とも共通骨格9 cluster、反応中心4 clusterである。
- ただしSBの低エネルギー領域には複数の近接構造があり、単一current-bestだけで比較するよりensembleで見る必要がある。電子E重みは振動・回転エントロピーを含まないため、存在比の予測値ではない。

## 仮説5：分子内相互作用／折り畳み

| 指標 | Bz | SB | SB−Bz |
| --- | ---: | ---: | ---: |
| benzoyl–core contact count | 2.071 ± 0.258 | 2.067 ± 0.250 | -0.004 |
| benzoyl–core contact score / Å | 1.556 ± 0.103 | 1.553 ± 0.099 | -0.003 |
| benzoyl ring–core contact count | 1.143 ± 0.515 | 1.134 ± 0.501 | -0.008 |
| TMS–core contact count | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 |
| nearest aromatic centroid distance / Å | 4.801 ± 0.017 | 4.808 ± 0.018 | 0.007 |
| heavy-atom radius of gyration / Å | 4.901 ± 0.053 | 5.632 ± 0.049 | 0.731 |
| carbonyl best clearance / Å | 0.286 ± 0.011 | 0.285 ± 0.010 | -0.001 |
| benzoyl–carbonyl torsion / deg | 40.141 ± 0.735 | 40.536 ± 0.659 | 0.395 |

- 共通骨格RMSD ≤0.50 Åの相互最近傍pairは4組。pair差の符号が揃わない指標はTMS固有効果と解釈しない。
- contactはvdW距離に基づく幾何学proxyで、安定化エネルギーではない。TMS接触数はBzに存在しない追加原子数の影響を受けるため、benzoyl ring–core共通部分と分けて扱う。
- N/O–H donorはなく、記録したC–H···O候補は弱い幾何学候補にすぎない。π配置もcentroid・面角の記述であり、π–π相互作用エネルギーを意味しない。

## 判定上の停止線

- 仮説4・5は『TMSが自由分子の配座分布や接触幾何を変え得るか』までを扱う。酵素内populationや反応速度には直接変換しない。
- Freqが部分集合のみなので、ここでの電子E重みを完全なBoltzmann populationとは呼ばない。
