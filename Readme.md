# Cross-species comparison of scRNA-seq

V0.0.1 (05/12/2023)  
Wrote by Yuichiro Hara

<br><br>


## 0. はじめに

複数の生物種のscRNA-seqにおいて、one-to-one orthologを用いて種間で遺伝子セットを共通させて統合解析を行う。

本稿では、
- Ensemblからのone-to-one orthologの取得
- 10x ChromiumによるscRNA-seqの実効
- Seuratによるデータ解析
  
という汎用性が高く簡便な方法を想定する。ただし、Ensemblにデータが格納されていない生物種では自身でオーソログ関係を推定する必要がある(別途記載予定)。
使用する生物をマウス、ニワトリ、スッポンとする。

<br><br>


## 1. 準備

### 1-1. 解析環境
#### R
Seurat (https://satijalab.org/seurat/articles/install_v5.html)  
Seuratで実行する統合解析用のライブラリも適宜インストールする  
Seurat-utils (Seurat <v5のみ; https://vertesy.github.io/Seurat.utils/)  

<br>

#### Python
Pandas  

<br>

Cell Ranger (by 10x Genomics; https://www.10xgenomics.com/jp/support/software/cell-ranger)
  
<br>

### 1-2. Cell Rangerの実効
Ensembl のオーソログアノテーションを用いる場合には、いずれの生物種も”同じリリースのEnsemblデータベース”からアセンブリと遺伝子アノテーションを取得し、リファレンスを作成する必要がある

#### 10x用リファレンスデータの作成
Ensembl(https://asia.ensembl.org/)からゲノムアセンブリ(FASTA)ファイルと遺伝子アノテーションファイル(GTF)を取得する。  
ページ最上部のDownloadに進み左側に現れるFTP downloadをクリック。https://asia.ensembl.org/info/data/ftp/index.html から直接アクセスできる。  
以前のリリースを用いたいときにはArchiveページ(https://asia.ensembl.org/info/website/archives/index.html)よりアクセスする。  
Celll Ranger mkrefコマンドを用いてリファレンスデータを作成する  
```
cellranger mkref --genome=PelSin_1.0.109 --fasta=Pelodiscus_sinensis.PelSin_1.0.dna_sm.toplevel.fa --genes=sP0.PelSin1ens109.ensinhouse.merged.ovlp.cat.gtf --memgb=24 --nthreads=16
```
上記ではEnsemblで配布された遺伝子アノテーションをそのまま用いているが、場合に応じてbulk RNA-seqなどを行いUTRの再アノテーションを行うこと。10x Chromiumの3' GEX kitを用いるならば、ChromiumのシーケンスリードでUTRを伸長させるPeak2UTRを用いるのも効果的。

全ての生物のscRNA-seqデータについてCell Rangerを実行しておく。例えば以下のような感じ  
```
cellranger count --id=turtle_cortex_1 --transcriptome=PelSin_1.0.109 --fastqs=/PATH/TO/fastq/ --memgb=24 --nthreads=16
``` 

<br><br>

## 2. One-to-one orthologの収集
EnsemblのBioMartにアクセス。Ensemblは定期的に新しいバージョンがリリースされるので、Cell Rangerに用いているゲノムにおけるEnsemblのリリースと一致しているか確認する。

用いる生物種の中で代表の生物種を決める。One-to-oneオーソログは種間で対称であるが、使用する遺伝子のシンボルを１つの生物種の命名に揃える必要があるため、その生物種を決めておく必要がある。
ここでは、マウス、ニワトリ、スッポンでone-to-oneを同定することを想定し、マウスを代表とさせる

BioMartで代表生物種を選択する。
“CHOSE DATABASE”でEnsemble Genes XXX (XXXはリリース番号)、
“CHOSE DATASER”で代表生物種(ここではMouse genes (GRCm39))を選択する。

ブラウザ上で左側のフレームから”Attributes”をクリックし、右フレームに現れた” Homologues (Max select 6 orthologues)”を選択する
“GENE”のカテゴリを開き、”Gene stable ID”, “Gene Name”の２つを選択する。  
“ORTHOLOGUES”のカテゴリ(アルファベット順に分かれている)を開き、対象とする生物種(ここでは”Chicken Orthologues”と” Chinese soft shell turtle Orthologues”)それぞれにおいて、”gene stable ID”, “gene name”, “homolog type”の３つを選択する。  

上部フレームのResultsをクリックし、TSVファイル(あるいはそのGZIP圧縮ファイル)として保存する。  
ここでは、`ortholog.mm-gg-ps.ens109.txt`としておく。以下の表がその抜粋。ファイルをexamplesフォルダに入れている。  

|Mouse gene stable ID|Mouse gene name|Chicken gene stable ID|Chicken gene name|Chicken homology type|Chinese softshell turtle gene stable ID|Chinese softshell turtle gene name|Chinese softshell turtle homology type|
|:---|:---|:---|:---|:---|:---|:---|:---
|ENSMUSG00000047161|Chst9|ENSGALG00010005573|CHST9|ortholog_one2one|ENSPSIG00000002647|CHST9|ortholog_one2one|
|ENSMUSG00000024304|Cdh2|ENSGALG00010005621|CDH2|ortholog_one2one|ENSPSIG00000002658|CDH2|ortholog_one2one|
|ENSMUSG00000087170|1700001G01Rik|||||||
|ENSMUSG00000048799|Cep120|ENSGALG00010001306|CEP120|ortholog_one2one|ENSPSIG00000014102|CEP120|ortholog_one2one|
|ENSMUSG00000073563|Csnk1g3||||ENSPSIG00000014585|CSNK1G3|ortholog_one2one|
|ENSMUSG00000059898|Dsc3|ENSGALG00010004256||ortholog_one2many|ENSPSIG00000002716||ortholog_many2many|
|ENSMUSG00000059898|Dsc3|ENSGALG00010004256||ortholog_one2many|ENSPSIG00000017548||ortholog_many2many|

<br>

スクリプト`get_onetoone.py`を用いてone-to-oneオーソログを抽出する。
```
python get_onetoone.py --sanitize ortholog.mm-gg-ps.ens98.txt > ortholog.mm-gg-ps.ens109.1to1.txt
```
<br>

出力は以下の形式のタブ区切りファイル。
|Mouse gene ID|Mouse gene name|Chicken gene ID|Chicken gene name|Chinese softshell turtle gene ID|Chinese softshell turtle gene name|
|:---|:---|:---|:---|:---|:---|
|ENSMUSG00000028101|Pias3|ENSGALG00000034042|PIAS3|ENSPSIG00000005069|PIAS3|
|ENSMUSG00000042354|Gnl3|ENSGALG00000001613|GNL3|ENSPSIG00000011486|GNL3|
|ENSMUSG00000027379|Bub1|ENSGALG00000008233|BUB1|ENSPSIG00000017705|BUB1|
|ENSMUSG00000035683|Melk|ENSGALG00000016438|MELK|ENSPSIG00000007833|MELK|
|ENSMUSG00000024346|Pfdn1|ENSGALG00000000946|PFDN1|ENSPSIG00000014687|PFDN1|

<br>


one-to-oneオーソログアノテーションの精度の問題で、着目したいone-to-oneオーソログがファイルに含まれないこともある。その場合には`--manualcuration`オプションでデータに含めるオーソログが記入されたファイルを指定する。
```
python get_onetoone.py --sanitize --manualcuration ortholog.whitelist.txt ortholog.mm-gg-ps.ens98.txt > ortholog.mm-gg-ps.ens109.1to1.curated.txt
```
このオプションで指定するファイルの形式は出力ファイルと同じ
|Mouse gene ID|Mouse gene name|Chicken gene ID|Chicken gene name|Chinese softshell turtle gene ID|Chinese softshell turtle gene name|
|:---|:---|:---|:---|:---|:---|
|ENSMUSG00000020950|Foxg1|ENSGALG00010006854|FOXG1|ENSPSIG00000000017|FOXG1|


  
10Xが提供するリファレンスゲノム(ヒト、マウス)を用いたいならば、他の生物種のリファレンスデータはEnsembl release 98から取得し、one-to-oneオーソログもrelease 98のBioMartより取得する必要がある。  

<br>

### OrthoFinderを用いた one-to-one orthologテーブルの作成
(執筆予定)

<br><br>

## 3. One-to-one orthologに絞ったscRNA-seqデータセットの作成
スクリプト`mkseuratobj.orthoset.prim.R`で、Cell Rangerの出力ファイルから、one-to-oneオーソログに遺伝子を絞り、かつ遺伝子名を代表の種(本稿ではマウス)のgene symbolに変換したSeuratオブジェクトを作成する。

```
Rscript mkseuratobj.orthoset.prim.R [CellRanger_output_dir] [One-to-one_table] [prefix] [species]
```
各引数は以下
```
[CellRanger_output_dir]: cellranger countで指定した出力のディレクトリ
[One-to-one_table]: 上記で作成したone-to-oneオーソログの対応表
[prefix]: 出力ファイルのprefix
[species]: 生物種名。one-to-oneオーソログの対応表の種名と同じにすること。スペースが含まれる場合にはダブルクオーテーションでくくる
```
出力は以下
```
[prefix].prim.rds CellRangerの出力をSeuratオブジェクト化したファイル。遺伝子数、遺伝子数ともそのまま
[prefix].ortho1to1.rds one-to-oneオーソログに遺伝子を絞り、かつ遺伝子名を代表の種のgene symbolに変換したSeuratオブジェクト
```
コマンド実行例はこちら
```
Rscript mkseuratobj.orthoset.prim.R /PATH/TO/turtle_cortex_1 ortholog.mm-gg-ps.ens109.1to1.txt turtle_cortex_1.1to1 "Chinese softshell turtle"
```

`mkseuratobj.orthoset.prim.R`はSeurat v3のAssayフォーマットに準拠しているため、v5のAssayフォーマットのオブジェクトを作成したい場合には`mkseuratobj.orthoset.prim.v5.R`を用いる

```
Rscript mkseuratobj.orthoset.prim.R /PATH/TO/turtle_cortex_1 ortholog.mm-gg-ps.ens109.1to1.txt turtle_cortex_1.1to1.v5 "Chinese softshell turtle"
```

作成したオブジェクトは種間で遺伝子数、遺伝子名が共通するため、この後は一般的な統合解析と同様に進められる。


<!--
4.	統合解析
-->