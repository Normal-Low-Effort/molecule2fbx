# 1Bz-LSD_RR / 1SB-LSD_RR 予備比較レポート

## 結論の範囲

この比較はB3LYP/def2-SVP・気相で得た計算候補構造の比較であり、実験構造、加水分解速度、酵素反応障壁を直接示すものではない。
現時点では、差が配座内分布または計算法の不確かさ以下なら『この条件では検出できない』とし、『効果が存在しない』とは扱わない。

- Bz: DFT代表 9 / 共通骨格unique 9、元のFreq選択 3/3完了、補足Freq 1、利用可能合計 4
- SB感度解析集合: DFT代表 11（元の10 + targeted follow-up 1） / 共通骨格unique 8、元のFreq選択 3/3完了、補足Freq 0、利用可能合計 6
- SB力場候補: 全重原子22 cluster → 共通骨格 9 cluster、DFT未カバー 1 cluster
- MBIS/CHELPG property対象: Bz 4、SB 6
- 過去のlauncher/startup失敗はoperational_historyとして保存し、最終Opt/Freq失敗とは別集計

主表は共通骨格clusterごとの最低電子エネルギー代表を1構造だけ採り、その集合を電子エネルギーで重み付けする。TMS回転だけの重複計上は避けるが、回転子の縮退や配置エントロピーを厳密に扱うものではなく、完全な分子ensembleのBoltzmann分布ではない。Freq未完了構造が残るため、Gibbs重みはFreq完了部分集合に対する条件付き値である。

## conformer別結果 — Bz

| conformer | ΔE / kJ mol⁻¹ | ΔG / kJ mol⁻¹ | Freq | 虚振動(<−20) | Opt由来 |
| --- | ---: | ---: | --- | ---: | --- |
| conf001 | 0.000 | 0.000 | 完了 | 0 | computed_original_ensemble_run |
| conf004 | 0.850 | 1.132 | 完了 | 0 | computed_original_ensemble_run |
| conf003 | 0.988 | 0.954 | 完了 | 0 | computed_original_ensemble_run |
| conf008 | 4.785 | 5.682 | 完了 | 0 | computed_original_ensemble_run |
| conf002 | 7.338 | — | 未実施 | null | computed_original_ensemble_run |
| conf007 | 7.630 | — | 未実施 | null | computed_original_ensemble_run |
| conf006 | 7.680 | — | 未実施 | null | computed_original_ensemble_run |
| conf005 | 7.786 | — | 未実施 | null | computed_original_ensemble_run |
| conf009 | 11.955 | — | 未実施 | null | computed_original_ensemble_run |

## conformer別結果 — SB

| conformer | ΔE / kJ mol⁻¹ | ΔG / kJ mol⁻¹ | Freq | 虚振動(<−20) | Opt由来 |
| --- | ---: | ---: | --- | ---: | --- |
| conf002 | 0.000 | 0.000 | 完了 | 0 | reused_external_read_only |
| conf007 | 1.006 | 0.896 | 完了 | 0 | computed_original_ensemble_run |
| conf005 | 1.015 | 1.156 | 完了 | 0 | computed_original_ensemble_run |
| conf004 | 1.091 | 0.603 | 完了 | 0 | reused_external_read_only |
| conf006 | 1.095 | 0.642 | 完了 | 0 | computed_original_ensemble_run |
| conf011_pool109 | 4.027 | — | 未実施 | null | computed_targeted_followup |
| conf003 | 7.187 | 6.618 | 完了 | 0 | reused_external_read_only |
| conf011 | 7.776 | — | 未実施 | null | computed_strict_selection_repair |
| conf010 | 7.824 | — | 未実施 | null | computed_original_ensemble_run |
| conf009 | 7.865 | — | 未実施 | null | computed_original_ensemble_run |
| conf008 | 7.919 | — | 未実施 | null | computed_original_ensemble_run |

## Bz–SB共通骨格の相互最近傍

| Bz | SB | 共通骨格RMSD / Å | 反応中心RMSD / Å |
| --- | --- | ---: | ---: |
| conf004 | conf004 | 0.742 | 0.810 |
| conf003 | conf005 | 0.042 | 0.010 |
| conf008 | conf002 | 0.496 | 0.003 |
| conf002 | conf009 | 0.773 | 0.006 |
| conf007 | conf011 | 0.015 | 0.005 |
| conf006 | conf010 | 0.036 | 0.012 |
| conf005 | conf008 | 0.012 | 0.005 |
| conf009 | conf003 | 0.505 | 0.006 |

原子対応はTMS末端枝を除去したグラフ同型写像で固定した。対称原子の全置換によるRMSD最小化は行っていないため、値は再現可能だが対称性を許した最小RMSDの上限になり得る。

## 指標比較

値は電子エネルギー重み付き平均 ± 配座間標準偏差。coverageは当該指標を持つ構造が全電子重みの何割を占めるかをBz / SBで示す。

| 指標 | Bz | SB | SB−Bz | weight coverage Bz/SB |
| --- | ---: | ---: | ---: | ---: |
| C=O結合長 / Å | 1.21559 ± 0.00026 | 1.21568 ± 0.00023 | 0.00010 | 1.000 / 1.000 |
| amide C–N結合長 / Å | 1.40226 ± 0.00194 | 1.40263 ± 0.00190 | 0.00037 | 1.000 / 1.000 |
| benzoyl–carbonyl二面角 / deg | 40.141 ± 0.735 | 40.536 ± 0.659 | 0.395 | 1.000 / 1.000 |
| 攻撃円錐accessible fraction | 0.1528 ± 0.0000 | 0.1528 ± 0.0000 | 0.0000 | 1.000 / 1.000 |
| 攻撃円錐best clearance / Å | 0.2860 ± 0.0108 | 0.2851 ± 0.0096 | -0.0009 | 1.000 / 1.000 |
| 攻撃円錐clearance p90 / Å | 0.0880 ± 0.0012 | 0.0883 ± 0.0012 | 0.0002 | 1.000 / 1.000 |
| MBIS carbonyl-C | 0.65068 ± 0.00042 | 0.64697 ± 0.00039 | -0.00371 | 0.929 / 0.933 |
| MBIS carbonyl-O | -0.51157 ± 0.00021 | -0.51134 ± 0.00021 | 0.00023 | 0.929 / 0.933 |
| CHELPG carbonyl-C | 0.53219 ± 0.01534 | 0.52892 ± 0.01322 | -0.00327 | 0.929 / 0.933 |
| CHELPG carbonyl-O | -0.48241 ± 0.00341 | -0.48465 ± 0.00230 | -0.00224 | 0.929 / 0.933 |
| Mayer C=O | 2.08854 ± 0.00466 | 2.08821 ± 0.00429 | -0.00032 | 1.000 / 1.000 |
| Mayer C–N | 1.03927 ± 0.00440 | 1.03792 ± 0.00419 | -0.00135 | 1.000 / 1.000 |
| dipole / D | 5.2785 ± 1.0022 | 5.3726 ± 1.1011 | 0.0941 | 1.000 / 1.000 |
| HOMO Loewdin population (benzoyl center) | 0.09320 ± 0.00190 | 0.09389 ± 0.00171 | 0.00069 | 0.929 / 0.933 |
| LUMO Loewdin population (benzoyl center) | 0.46325 ± 0.00034 | 0.44171 ± 0.00067 | -0.02154 | 0.929 / 0.933 |
| HOMO energy / eV | -5.5334 ± 0.0438 | -5.5115 ± 0.0429 | 0.0218 | 1.000 / 1.000 |
| LUMO energy / eV | -1.4603 ± 0.0339 | -1.4940 ± 0.0337 | -0.0337 | 1.000 / 1.000 |
| carbonyl stretch / cm⁻¹ | 1777.46 ± 0.28 | 1777.01 ± 0.36 | -0.45 | 0.929 / 0.880 |

## 仮説1（Si–C結合・配置による立体障害緩和）

### 支持と整合する点

- 固定Bürgi–Dunitz円錐でのaccessible fractionはBz/SBとも0.1528、best clearanceのSB−Bzは−0.0009 Åで、配座内標準偏差（約0.01 Å）よりはるかに小さい。少なくとも静的気相構造では、SBのpara-TMSがcarbonyl Cへの幾何学的アクセスを大きく閉じる差は検出されなかった。
- SBのcarbonyl C–Si距離は電子エネルギー重み付きで約6.25 Åで、TMS重原子はcarbonyl Cから5 Å以内に0個だった。TMSは見かけの体積ほど反応中心を直接覆っていないという考えとは整合する。
- 共通骨格RMSDにより、TMSメチル回転だけを別のLSD–benzoyl配座として数える問題を除いた。

### 反する可能性のある点

- SBでclearanceが増えた証拠はなく、平均はごくわずかに低い。したがって『Si–C結合が長いこと自体がアクセスを改善する』までは支持できない。現状が示すのは、para-TMSによる直接遮蔽をこの指標では検出できない、という範囲である。
- accessible fractionは全構造で同一になり、角度刻みと閾値に量子化されている。best clearance、p90、probe-radius感度を併記しても、これは予備的な幾何学指標である。
- この幾何学probeには水、酵素ポケット、基質誘導適合、遷移状態が含まれない。静的アクセスが同じでも加水分解障壁が同じとは限らない。

## 仮説2（TMSの電子効果）

### 支持と整合する点

- carbonyl-CのMBIS電荷はSB−Bz = −0.00371 eで、各分子内の配座標準偏差（約0.0004 e）より大きく、低エネルギーproperty部分集合では一貫した差として検出された。TMS-Bz側でcarbonyl Cがわずかに低正電荷になる方向である。
- benzoyl centerのLUMO Loewdin populationはSB−Bz = −0.02154で、配座標準偏差より大きい。C=Oが+0.00010 Å長く、Mayer C=Oが−0.00032、carbonyl stretchが−0.45 cm⁻¹という小さな変化も、弱い電子供与を想定した方向とは概ね整合する。

### 反する可能性のある点

- CHELPG carbonyl-C差（−0.00327 e）は配座標準偏差（約0.013–0.015 e）より小さく、結合長、Mayer bond order、dipole、stretchの差も多くは配座ばらつき以下である。MBIS/LUMO局在化ほど全指標が強く一致しているわけではない。
- したがって『TMSに電子効果がある』という予備的証拠は得たが、その大きさが溶液中でも維持されるか、加水分解を速めるか遅くするかは判断できない。電子差の検出と反応速度の説明は分ける。

## SB未DFTクラスタの再評価

| pool index | 安価な電子構造screen ΔE / kJ mol⁻¹ | 判定 |
| ---: | ---: | --- |
| 46 | 6.481 | review_for_one_additional_opt |

pool 109は事前記録後に1構造だけtargeted DFT Optし、既存current best conf002より4.027 kJ mol⁻¹高かった。共通骨格RMSD 0.75 Åでは既存構造と重複せず、RRを保持したが、Freq未実施のためimaginary_modesはnullである。
残るpool 46はMMFF構造上のDFT単一点で+6.481 kJ mol⁻¹、GFN2-xTB Optで+6.882 kJ mol⁻¹だった。未緩和単一点/別モデルの順位なので、追加Opt候補の優先順位付けにだけ用いる。今回は大量投入せず保留した。

## 現時点で判断できないこと

- Bz→SBの置換では立体効果と電子効果が同時に変わるため、観測差の因果分離はこの二分子だけでは完結しない。
- Freqで虚振動が0でも、それは当該構造が局所極小候補であることを示すだけでglobal minimumを保証しない。
- 低周波モードを含む調和近似のΔGはsub-kJ/mol順位に敏感であり、Freq部分集合だけの平均を完全ensemble平均とは呼べない。

## 探索範囲の非対称性

Bzは最大10枠に達する前に候補が収まった一方、SBは元の全重原子cluster数が22で最大10枠に切られた。共通骨格では9 clusterに縮約し、pool 109を追加Optした結果、未カバーはpool 46の1 clusterまで減った。ただしBz/SBの初期候補選択履歴は完全対称ではないため、SB側current bestの確度はなおわずかに低い。

## 次に価値が高い計算

1. 直ちに追加の長時間DFTは不要。まず固定円錐より連続的な求核攻撃trajectory/SASA指標を実装し、相互最近傍のBz/SB配座でpaired差を確認する（新規量子化学計算なし）。
2. 仮説2のMBIS/LUMO差を検証するなら、相互対応する低エネルギー構造だけに高い基底・dispersion・暗黙溶媒を用いたproperty単一点を行う。OptやTSより安価で、気相/基底依存性を直接判定できる。
3. pool 46 Optは探索非対称性を完全に閉じたい場合だけ実行する。現screenではcurrent bestを更新する優先度は低く、見積り0.4–0.8時間、上位化した場合のみFreq約2時間を追加する。
4. 上記で電子差または立体差がモデル変更後も維持されて初めて、溶媒中反応物複合体、加水分解TS、置換基対照へ進む価値を再判定する。
