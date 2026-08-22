# molecule2fbx 0.4.0

[英語版](README.md)

`molecule2fbx` は、PubChem CIDまたはSMILESから球棒分子モデルを生成し、
Blenderを介してFBXへ出力するコマンドラインツールです。PubChem 3D構造、
RDKit力場構造、ORCA量子化学構造を区別し、由来と計算条件をJSONメタデータへ
保存します。
このリポジトリには、molecule2fbxを用いて実施された研究解析の厳選されたスナップショットも含まれています。大規模な生計算データはローカルに保存されており、リポジトリには含まれていません。研究データにご興味がある場合は、お気軽にお問い合わせください。

バージョン0.4.0では、ETKDGによる配座生成、DFT構造最適化、DFT後の重複除去、
選択的な振動数計算、相対エネルギーとGibbs自由エネルギーの解析を一つの処理として
実行する `--ensemble` モードを追加しました。

計算構造は実験構造ではありません。一回の計算で得られた最低エネルギー構造も、
調べた候補中の現時点での最良構造にすぎず、真の大域的最小構造を保証しません。

## 既存アンサンブルの比較解析

`molecule2fbx.comparison` は、従来の全重原子クラスタリングを置き換えず、
追加の解析指標を提供します。

- `all_heavy_rmsd_cluster_id`：入力時の原子順を固定した全重原子RMSD。
- `common_scaffold_rmsd_cluster_id`：末端の芳香環-Si(CH3)3分岐、すなわちSiと
  三つのメチル炭素を除外し、入力時の原子順を固定したRMSD。
- `reaction_center_rmsd_cluster_id`：芳香族N-ベンゾイルカルボニル炭素から
  グラフ上で二結合以内にある重原子を対象としたRMSD。

対称原子の入れ替えは行いません。ORCAのXYZは入力時の原子順を保持する一方、
自動的な原子置換は原子の同一性を暗黙に変える可能性があるためです。共通骨格指標では
末端TMS分岐全体を除き、等価なメチル基の回転や入れ替えだけで別のLSD-ベンゾイル
骨格クラスタが生じないようにしています。この指標は追加の記述子であり、座標や従来の
全重原子クラスタIDを書き換えません。

Bz体とTMS-Bz体を直接比較するRMSDでは、RDKitグラフから末端TMS分岐を除き、
キラリティを保持するグラフ同型写像を求めます。ベンゾイル基のC/O/N/ipso原子を
固定点とし、保持された入力原子順で変位が最小となる写像を選び、その一つの決定論的な
写像にKabsch重ね合わせを適用します。RMSDを下げるためのベンゼン環やエチル分岐の
対称置換探索は行いません。採用した原子対応とこの仮定は `analysis.json` に記録します。
したがって分子間RMSDは再現可能で原子同一性を保ちますが、未探索の対称写像によって
さらに低いRMSDが得られる場合には上限値となり得ます。

記述子レポートには、従来の全代表構造による重み付けに加え、各共通骨格クラスタから
電子エネルギー最低の構造を一つだけ数える `common_scaffold_unique` を収録します。
Bz/SB感度比較では、TMSだけの回転重複を独立した骨格配座として数えないため、後者を
主に使用します。これは回転子の縮退や配置エントロピーを厳密に扱う方法ではなく、
その制限を明記しています。

Bz/SB予備解析は次のコマンドで再生成できます。

```powershell
work\test-venv\Scripts\python.exe scripts\analyze_bz_sb_preliminary.py
```

この解析は、振動数計算がない場合に `imaginary_modes: null` とし、最終的な計算失敗と
再開時の回復履歴を分離します。また、Opt/Freqの由来を記録し、決定論的な力場配座群を
再生成し、保存されたORCA物性を抽出し、サンプリングしたBurgi-Dunitz円錐上で方向別の
カルボニル接近性を評価します。立体プローブは反応物構造の記述子にすぎず、酵素、溶媒、
活性化障壁を表すモデルではありません。最終DFT代表構造のすべてにFreqがない場合、
Gibbs重み付き値は条件付きと明記します。

## リポジトリと研究データ

Gitで管理するファイルと、ローカルに保持する大容量計算データを分離しています。

```text
molecule2fbx/      コマンドラインツール本体
tests/             自動テスト
scripts/           アンサンブル実行・追加計算・解析スクリプト
research_results/  Git共有用の軽量な解析結果スナップショット
outputs/           Git対象外のORCA・xTB・runEDDB・FBX記録
work/              仮想環境・試験・パッケージ作業領域
tools/             Git対象外のローカル外部ツール
```

`outputs/` が計算記録の正本です。`research_results/` は
`scripts/export_research_results.py` が明示的に選択してコピーした小容量の要約、表、
現時点での最良構造のXYZです。ORCA原データの代替ではありません。

```powershell
python scripts\export_research_results.py --check
python scripts\export_research_results.py
```

スナップショット中のローカルなワークスペース絶対パスは `${WORKSPACE}` へ置換されます。
各ファイルのコピー元、サイズ、原本とスナップショット双方のSHA-256は
`research_results/SNAPSHOT_MANIFEST.json` に記録されます。

## 必要環境

- Python 3.9以上
- Blender 3.xまたは4.x
- RDKit、NumPy、requests
- DFTまたはHF計算を行う場合は、別途インストールしたORCA 6.x

```powershell
python -m pip install .
```

開発環境：

```powershell
python -m pip install -e ".[dev]"
python -m pytest
```

ORCAはライセンスと配布条件のため同梱していません。`--orca`、
`ORCA_EXECUTABLE`、`ORCADIR`、または `PATH` で指定します。

```powershell
$env:ORCA_EXECUTABLE = "C:\Orca_6.1.1\orca.exe"
$env:BLENDER_EXECUTABLE = "C:\Program Files\Blender Foundation\Blender 4.0\blender.exe"
```

## 既存のコマンドライン操作

従来の位置引数によるCID入力はそのまま使用できます。

```powershell
molecule2fbx 5360697
molecule2fbx --cid 5360697 --method auto
```

`auto` はPubChem 3Dを優先し、利用できない場合だけETKDGとMMFF94sまたはUFFへ
切り替えます。`auto` が量子化学計算を開始することはありません。

```powershell
molecule2fbx --smiles "CCO" --name Ethanol --method forcefield
molecule2fbx --smiles "O" --name Water --method dft
```

ORCAの既定電子構造条件は `B3LYP/def2-SVP`、電荷0、多重度1です。

```powershell
molecule2fbx --smiles "[O]" --method dft `
  --functional B3LYP --basis "6-31G(d)" `
  --charge 0 --multiplicity 3
```

`--nprocs` を省略すると `os.cpu_count()` を使用します。WindowsがRyzen 5 5600を
12論理プロセッサとして認識していれば、ORCA入力には `%pal nprocs 12 end` が生成されます。

## 配座アンサンブルモード

完全な立体化学を含むSMILESを指定します。`--ensemble` 自体が長時間計算への明示的な
許可となるため、このモードで `--method` を省略するとDFTを選択します。

```powershell
molecule2fbx --ensemble `
  --smiles "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C(=O)c4ccc(cc4)[Si](C)(C)C)c5cccc(C2=C1)c35" `
  --name 1SB-LSD_RR `
  --nprocs 16 --maxcore 1000
```

既定の処理：

```text
SMILES・電荷・多重度・立体化学を検証
  → ETKDGで200構造を生成
  → MMFF94s、使用不能ならUFFで最適化
  → 力場最低エネルギーから10 kJ/mol以内を保持
  → 重原子RMSD 0.75オングストロームでクラスタリング
  → 最大10代表をORCA B3LYP/def2-SVPで最適化
  → DFT構造を0.75オングストロームで再クラスタリング
  → 電子エネルギー最低構造から5 kJ/mol以内の最大3構造でFreq
  → 振動数・ZPE・熱補正・H・G・相対E・相対Gを保存
```

主要なオプション：

| オプション | 既定値 | 内容 |
|---|---:|---|
| `--conformer-pool` | 200 | ETKDG生成数。アンサンブルでは100～500 |
| `--conformers` | 10 | DFTへ送るクラスタ代表の上限 |
| `--random-seed` | 61453 | ETKDG乱数シード |
| `--embedding-prune-rmsd` | -1 | 埋め込み時の枝刈り。-1は無効 |
| `--forcefield-energy-window-kj` | 10 | DFT前の力場エネルギー窓 |
| `--conformer-rmsd-threshold` | 0.75オングストローム | DFT前の重原子RMSD閾値 |
| `--dft-energy-window-kj` | 制限なし | 任意のDFT後エネルギー窓 |
| `--dft-rmsd-threshold` | 0.75オングストローム | DFT後の重原子RMSD閾値 |
| `--frequency-window-kj` | 5 | Freq選択の電子エネルギー窓 |
| `--frequency-max` | 3 | 自動選択するFreq構造の上限 |
| `--frequency-include N` | なし | 高エネルギーでも配座NをFreqへ保持。反復可 |
| `--imaginary-threshold-cm1` | -20 | これ未満を虚振動として判定 |
| `--low-frequency-threshold-cm1` | 50 | これ未満を低周波として記録 |

RMSD閾値を自動的に緩和しません。独立クラスタが10個未満なら、得られた代表だけを
DFTへ送ります。

RDKitのMMFF/UFF一括最適化機能は、Python APIから実際の反復回数を返しません。
そのためメタデータには `forcefield_optimization_iterations: null` と
`forcefield_iteration_count_available: false` を保存し、状態と反復上限を併記します。
反復回数を推測して記録することはありません。

## 立体化学

通常モードでは未指定立体中心を警告します。`--strict-stereochemistry` または
`--ensemble` では、未指定立体中心やE/Z候補があればORCA開始前に停止します。

配座探索と立体異性体探索は別処理です。入力SMILESの `@` と `@@` を保持し、未指定中心から
R/S異性体を列挙しません。正準異性体SMILES、CIPラベル、未指定中心をメタデータへ保存します。

## 既存ORCA計算の再利用

```powershell
molecule2fbx --ensemble --smiles "..." --name Molecule `
  --reuse-calculations output\Molecule_dft_calculations `
  --nprocs 16 --maxcore 1000
```

再利用時は、ORCA入力・出力・XYZ、正常終了、Opt収束、原子順、正準異性体SMILES、
計算法、汎関数、基底関数、電荷、多重度を検証します。互換なOptは再実行せず、不足する
配座だけを計算します。不完全または条件不一致の既存ディレクトリは上書きしません。

現在のアンサンブル実装は気相最適化で、追加の分散補正キーワードを使用しません。
これらの条件はアンサンブルJSONへ記録し、異なる条件の結果を暗黙に一つのアンサンブルとして
比較しません。

既存の最適化済みXYZへFreqだけを追加する場合：

```powershell
molecule2fbx --frequency-only path\to\conformer_001.xyz `
  --nprocs 16 --maxcore 1000
```

この入力にはFreqを含みますが、Optと `%geom` は含みません。電荷、多重度、汎関数、
基底関数をメタデータまたは対応する `.inp` から復元します。Freq専用解析ではOpt収束印を
要求しません。

## 出力

アンサンブル出力例：

```text
output/
├─ Molecule_dft.fbx
├─ Molecule_dft.metadata.json
├─ Molecule_dft_conf001.fbx
├─ Molecule_dft_conf001.metadata.json
├─ Molecule_dft_ensemble.json
└─ Molecule_dft_calculations/
   ├─ conformer_001/
   └─ frequency_additions/
```

アンサンブルJSONには、力場候補の除外理由、初期クラスタ、全DFT収束構造の相対電子
エネルギー、DFT後クラスタと重複先、最終代表、Freq対象、全振動数、低周波、虚振動、
熱力学量、相対Gibbs自由エネルギー、再利用状態、計算ディレクトリ、ソフトウェア条件、
時刻を保存します。

Freq未実行構造を局所極小確認済みとは扱いません。判定は次の三種類です。

- `local_minimum_candidate`：Freq完了、-20 cm-1未満の虚振動なし。
- `not_a_confirmed_local_minimum`：虚振動あり。
- `not_evaluated`：Freq未実行または解析不能。

低周波モードは調和振動近似によるGibbs自由エネルギーへ大きく影響する可能性があります。
ORCAがQuasi-RRHOを使用した場合はメタデータへ記録しますが、熱力学量の精度や化学的な
正しさを保証するものではありません。

## 金属と特殊分子

金属を検出した場合、有機分子向け既定値を自動適用しません。計算には
`--allow-metals` と、汎関数、基底関数、電荷、多重度の明示が必要です。

```powershell
molecule2fbx --smiles "..." --method dft --allow-metals `
  --functional B3LYP --basis def2-TZVP --charge 0 --multiplicity 1
```

この保護機能は、基底関数、ECP、酸化状態、開殻状態が適切であることを保証しません。
金属錯体では専門的な電子状態の検討が必要です。

## 構成

- `config.py`：コマンドライン画面と将来のGUIで共有する要求・検証モデル。
- `cli.py`：コマンドライン引数の解析。
- `structures.py`：SMILES検証、ETKDG、MMFF94s/UFF、重原子RMSD。
- `ensemble.py`：エネルギー選別、クラスタリング、相対エネルギー、アンサンブルJSON。
- `pipeline.py`：利用者画面に依存しない処理統合。
- `quantum/base.py`：バックエンド非依存の設定、結果、熱力学データ。
- `quantum/orca.py`：ORCA入力生成、外部実行、出力解析。
- `quantum/reuse.py`：OptとFreqの非破壊的な再利用。
- `frequency.py`：最適化済みXYZに対するFreq専用処理。
- `blender_export.py` と `blender_worker.py`：FBX生成。

ORCAを採用した理由は、Windows上で外部実行ファイルとして導入でき、DFT構造最適化、
開殻系、振動数計算、広い基底関数・ECP選択を扱え、大きな量子化学実行環境をPython環境へ
埋め込まずに済むためです。

`outputs/1SB-LSD_RR_redo` はローカルの検証データとして保持します。そのconf002は既存の
四構造中で現時点の最良構造ですが、大域的最小構造ではありません。より広い配座探索で、
同等の条件においてさらに低いエネルギー構造が得られた場合のみ、新しい最良候補とします。

注:英語があまり得意ではないため、READMEの作成にはAIを活用しました。また、プロジェクト全体を通じて、開発の補助としてもAIを利用しています。ただし、研究の設計、計算、および最終的な判断は、私自身が行ったものです。
