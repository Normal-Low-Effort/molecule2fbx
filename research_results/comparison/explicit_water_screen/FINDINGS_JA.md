# 仮説6 明示的1水分子microsolvation screen

対応する2組のBz/SB配座について、neutral/N6H+それぞれのN1-carbonyl周辺へ水1分子を9方向から配置し、GFN2-xTB/ALPB(水)で最適化した。これは明示的溶媒MDでも加水分解反応計算でもない。O–Hが1.30 Åを超えた配置は水のproton transfer/dissociationとして局所水和集計から除外した。

| pair | state | Bz carbonyl-local/intact | SB carbonyl-local/intact | Bz best Eint | SB best Eint | SB−Bz / kJ mol⁻¹ |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| pair_low_energy | neutral | 1/7 | 1/7 | -7.98 | -9.61 | -1.63 |
| pair_low_energy | n6_protonated | 1/7 | 1/7 | -9.69 | -10.03 | -0.35 |
| pair_close_geometry | neutral | 1/7 | 2/7 | -9.42 | -9.72 | -0.29 |
| pair_close_geometry | n6_protonated | 1/7 | 1/7 | -15.25 | -15.48 | -0.24 |

Eintは最適化複合体座標で E(complex)−E(solute fragment)−E(water fragment) としたxTB相互作用energy proxy。溶質変形は別項に分離したが、標準状態、エントロピー、濃度、BSSE、多水分子効果を含まないため結合自由エネルギーとは呼ばない。

このscreenでcarbonyl-local保持率と対応pairの差が一貫しない場合、TMSによる局所第一水和殻効果はこの解像度では検出不能とする。差が一貫しても、十分な水和ensembleと反応障壁で再検証するまで速度差へ結び付けない。
