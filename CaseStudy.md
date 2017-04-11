[生命動態のデータサイエンス](https://github.com/haruosuz/DS4GD)

----------

# Case Study
**ケーススタディ**

目次
- [Unixコマンド]
- [UniProtKB Swiss-Prot protein sequence database](#uniprotkb-swiss-prot-protein-sequence-database)
- [Silva rRNA database](#silva-rrna-database)
- [NCBI assembly summary](#ncbi-assembly-summary)

----------

### Unixコマンド
![http://techacademy.jp/magazine/5155](http://static.techacademy.jp/magazine/wp-content/uploads/2015/01/ss-1-620x375.jpg)

[ターミナル](http://techacademy.jp/magazine/5155)を開き、`bash`を起動する:  

	bash

[ターミナルの主要コマンド](https://techacademy.jp/magazine/5155#sec3)

	# ディレクトリ間を移動する
	cd

	# 現在のディレクトリ内に存在する、ディレクトリやファイルを表示する
	ls

	# 今開いているディレクトリまでのパスを表示する
	pwd

	# 「projects」ディレクトリを作成する
	mkdir projects

	# 「projects」ディレクトリへ移動する
	cd projects/

	# 空ファイルを作成する
	touch test.txt

	# ファイルをコピーする
	cp test.txt test2.txt

	# ファイルやディレクトリの名前を変更する
	mv test.txt test1.txt

	# ファイルを削除する
	rm test2.txt

	# コマンドのマニュアルを開く。終了する場合に「q」を入力
	man man

	# ファイルやディレクトリの検索
	find .

	# ターミナル画面をクリアし、表示された文字を全て消去する
	clear

	# ログアウトする
	exit

----------

## UniProtKB Swiss-Prot protein sequence database
[Swiss-Prot](https://ja.wikipedia.org/wiki/Swiss-Prot)タンパク質配列データベース

### Website
[UniProt](http://www.uniprot.org)の[Download latest release](http://www.uniprot.org/downloads)を開く。
<ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/> をブラウザ（Firefox または Chrome）で開く。
*uniprot_sprot.fasta.gz* を右クリックし、「リンクのURLをコピー (Copy Link)」する。

### Download

![http://techacademy.jp/magazine/5155](http://static.techacademy.jp/magazine/wp-content/uploads/2015/01/ss-1-620x375.jpg)

[ターミナル](http://techacademy.jp/magazine/5155)を開く。

プロジェクト・ディレクトリを作成し移動する:  

    mkdir -p ~/projects/uniprot_sprot/data
    cd ~/projects/uniprot_sprot/data/

[FASTA](https://ja.wikipedia.org/wiki/FASTA)形式の圧縮ファイル（*uniprot_sprot.fasta.gz*）を`wget`または`curl`でダウンロードする:  

    wget ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz
    # or
    curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz

### Inspecting Data

`ls -l`でファイルの詳細情報を表示する:  

    ls -l

`gunzip`コマンドで解凍する:  

    gunzip -c uniprot_sprot.fasta.gz > uniprot_sprot.fasta

`ls -lh`でファイルサイズを確認する:  

    ls -lh

`head`コマンドを用いて、ファイルの先頭部分を表示する:  

    head uniprot_sprot.fasta

[パイプ](https://ja.wikipedia.org/wiki/パイプ_%28コンピュータ%29)でプログラムの入出力をつなぐ。

[FASTA headers](http://www.uniprot.org/help/fasta-headers)

`grep`コマンドを用いて、FASTAファイルのヘッダ（`>`で始まる行）にマッチする行を抽出し、`head`で先頭部分を表示する:  

    grep '^>' uniprot_sprot.fasta | head

配列の数をカウントする:  

    grep '^>' uniprot_sprot.fasta | wc -l

[アミノアシルtRNA合成酵素](https://ja.wikipedia.org/wiki/アミノアシルtRNA合成酵素) (aminoacyl-tRNA synthetase/synthase) を対象として、文字列'tRNA synth'にマッチする行を抽出し、カウントする:  

    grep 'tRNA synth' uniprot_sprot.fasta | wc -l

[リボソーム](https://ja.wikipedia.org/wiki/リボソーム)タンパク質 ('ribosomal protein') にマッチする行をカウントする:  

    grep -ic 'ribosomal protein' uniprot_sprot.fasta

`grep`コマンドは、`-c`オプションでパターンにマッチした行数を表示し、`-i`オプションで大文字小文字を区別しない（ignore case）。

[ベロ毒素](https://ja.wikipedia.org/wiki/ベロ毒素)([Shiga toxin](https://en.wikipedia.org/wiki/Shiga_toxin)) を検索する。'Shiga'にマッチする行を表示する:  

    grep '^>' uniprot_sprot.fasta | grep 'Shiga'

### Working with Data in R

[R の起動と終了](http://cse.naro.affrc.go.jp/takezawa/r-tips/r/02.html)  

![http://cse.naro.affrc.go.jp/takezawa/r-tips/r/02.html](http://cse.naro.affrc.go.jp/takezawa/r-tips/r/image/Mac.gif)

[作業ディレクトリ](http://cse.naro.affrc.go.jp/takezawa/r-tips/r/06.html)の変更と確認:  

    setwd('~/projects/uniprot_sprot/data/')
    getwd()
    dir()

[パッケージ](http://cse.naro.affrc.go.jp/takezawa/r-tips/r/08.html)`seqinr`を呼び出す:  

    library(seqinr)

`read.fasta()`関数で配列データを読み込む:  

    ld <- read.fasta(file = 'uniprot_sprot.fasta.gz', seqtype = c('AA'), strip.desc = TRUE)

配列の数をカウントする:  

    length(ld)

- [文字列 | Rを利用して文字列のマッチング,結合,分割,置換を行う関数](http://stat.biopapyrus.net/r/string.html)

`getAnnot`関数を用いて、配列のアノテーションを取得する:  

    annotation <- getAnnot(ld)

`grep(pattern, x)`は、`pattern`にマッチするベクトル`x`の全要素の番号を返す。
'Alpha-actinin-2'にマッチする要素番号を取得:  

    i <- grep(pattern='Alpha-actinin-2', x=annotation, ignore.case=TRUE)

[リスト](http://cse.naro.affrc.go.jp/takezawa/r-tips/r/23.html)の成分を取り出す:  

    ld[i]

`write.fasta()`関数を用いて、配列データをFASTA形式ファイルとして書き出す:  

    write.fasta(sequences=ld[i], names=getName(ld[i]), file.out='myseq.fasta')

作業を中断し再開する（Rを終了し再起動する）。作業ディレクトリを変更し、パッケージ`seqinr`を呼び出し、`read.fasta()`関数で配列データを読み込む:  

    setwd('~/projects/uniprot_sprot/data/')
    library(seqinr)
    lx <- read.fasta(file = 'myseq.fasta', seqtype = c('AA'), strip.desc = TRUE)
    unlist(getAnnot(lx))

`sapply()`は、リストの各要素に関数を適用する。アミノ酸残基数を求める:  

    sapply(lx, length)

#### [Pairwise Sequence Alignment](https://github.com/haruosuz/r4bioinfo/tree/master/R_Avril_Coghlan#pairwise-sequence-alignment)
ペアワイズ配列アラインメント

[ドットプロットで2つの配列を比較](https://github.com/haruosuz/r4bioinfo/tree/master/R_Avril_Coghlan#comparing-two-sequences-using-a-dotplot)

    s1 <- lx[[1]]
    s2 <- lx[[3]]
    dotPlot(s1, s2)
    wsize <- 3; dotPlot(s1, s2, wsize=wsize, wstep=wsize, nmatch=wsize, xlab=getName(s1), ylab=getName(s2))

for文でループ処理を行う。複数のドットプロットをファイル出力する:  

    pdf("Rplot.pdf") # create a PDF device called Rplot.pdf
    for (n1 in 1:length(lx)) {
      for (n2 in n1:length(lx)) {
        s1 <- lx[[n1]]
        s2 <- lx[[n2]]
        dotPlot(s1, s2, xlab=getName(s1), ylab=getName(s2))
      }
    }
    dev.off(which = dev.cur()) # close the device

[2つのタンパク質配列間のグローバル・アライメント](https://github.com/haruosuz/r4bioinfo/tree/master/R_Avril_Coghlan#pairwise-global-alignment-of-protein-sequences-using-the-needleman-wunsch-algorithm)

    # Biostringsパッケージの呼び出し
    library(Biostrings)
    # データ(置換行列)のロード
    data(BLOSUM50)

    # 文字ベクトルを文字列に変換（大文字に変換）し、pairwiseAlignment
    aln_global <- pairwiseAlignment(toupper(c2s(s1)), toupper(c2s(s2)), substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
    aln_global
    writePairwiseAlignments(aln_global)
    writePairwiseAlignments(aln_global, file="aln_global.txt")

[2つのタンパク質配列間のローカル・アライメント](https://github.com/haruosuz/r4bioinfo/tree/master/R_Avril_Coghlan#pairwise-local-alignment-of-protein-sequences-using-the-smith-waterman-algorithm)

    aln_local <- pairwiseAlignment(toupper(c2s(s1)), toupper(c2s(s2)), substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE, type="local")
    aln_local
    writePairwiseAlignments(aln_local, file="aln_local.txt")

#### [Multiple sequence alignment](https://github.com/haruosuz/r4bioinfo/tree/master/R_msa#multiple-sequence-alignment)
多重配列アライメント (多重整列) 

    library(Biostrings)
    aaseq <- readAAStringSet(file = 'myseq.fasta')

    library(msa)    aln <- msa(aaseq, method="Muscle")
    aln

    # 多重配列アライメントをファイル出力
    writeXStringSet(unmasked(aln), file="aln.fasta")

[SeaView](http://www2.tba.t-com.ne.jp/nakada/takashi/phylogeny/seaview2.html)を用いて、FASTA形式ファイルの多重配列アライメントを表示する。

![](http://www2.tba.t-com.ne.jp/nakada/takashi/phylogeny/fig/sv2/sv7.jpg)

[WebLogo](http://blog.amelieff.jp/?eid=210264)を用いて、配列の保存度を表示する。

![](http://blog.amelieff.jp/images/weblogo2.png)

#### [Phylogenetic analysis](https://github.com/haruosuz/r4bioinfo/tree/master/R_msa#phylogenetics)
系統解析 

    library(phangorn) # read.aa # dist.ml

    # アミノ酸配列（FASTA形式ファイル）を読み込む
    aln <- read.aa(file = "aln.fasta", format = "fasta")

    # 配列間の距離を計算
    d <- dist.ml(aln, model="Dayhoff")

    # 近隣結合法 NJ (Neighbor-Joining)
    tre <- nj(d)

    # 系統樹を描画
    plot.phylo(tre, type = "phylogram") # type "phylogram", "cladogram", "fan", "unrooted", "radial"

    # 系統樹をnewick形式でファイル出力
    write.tree(tre, file="myseq.newick")

[SeaView](http://doua.prabi.fr/software/seaview)または[FigTree](http://www.geocities.jp/ancientfishtree/FigTree.html)を用いて、newick形式ファイルの系統樹を表示する。

![](http://en.bio-soft.net/tree/figtree.jpg)

#### Amino acid usage
タンパク質のアミノ酸組成

- [Lobry, J.R., Chessel, D. (2003) Internal correspondence analysis of codon and amino-acid usage in thermophilic bacteria. Journal of Applied Genetics, 44:235-261.](http://pbil.univ-lyon1.fr/members/lobry/repro/jag03/)
- [Lobry JR, Gautier C. Nucleic Acids Res. 1994 Aug 11;22(15):3174-80. 'Hydrophobicity, expressivity and aromaticity are the major trends of amino-acid usage in 999 Escherichia coli chromosome-encoded genes.'](http://www.ncbi.nlm.nih.gov/pubmed/8065933)

`AAstat()`関数を用いて、タンパク質の配列情報（アミノ酸残基数、物理化学的クラスの割合、等電点の理論値）を求める:  

    AAstat(lx[[1]])

アミノ酸の個数を出力:  

    AAstat(lx[[1]], plot = FALSE)$Compo

物理化学的クラス (Tiny, Small, Aliphatic, Aromatic, Non-polar, Polar, Charged, Positive, Negative) のうち、[芳香族アミノ酸](https://ja.wikipedia.org/wiki/芳香族アミノ酸) Aromatic の割合を出力:  

    AAstat(lx[[1]], plot = FALSE)$Prop$Aromatic

`sapply()`関数は、リストの各要素に関数を適用する。複数タンパク質配列の物理化学的クラスの割合を求める:  

    sapply(lx, function(x) AAstat(x, plot=FALSE)$Prop )

複数タンパク質配列のアミノ酸使用の絶対度数と相対度数を求める:  

    # 絶対度数
    sapply(lx, function(x) AAstat(x, plot=FALSE)$Compo )

    # 相対度数
    X <- sapply(lx, function(x) summary(x)$composition )

データをカンマ区切りファイルとして出力する:  

    write.csv(t(X), file="table.csv")
タンパク質のアミノ酸使用パターンを比較する。`matplot()`関数でプロットする:  

    matplot(X, type="l", col=rainbow(12), lty=1:6)
    legend("topleft", legend=colnames(X), col=rainbow(12), lty=1:6)
    # lty: 0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash

[クラスター分析](https://github.com/haruosuz/DS4GD/blob/master/hclust.md)

    plot(hclust(dist(t(X))))

[ヒートマップ](https://github.com/haruosuz/DS4GD/blob/master/hclust.md#heat-map)

    heatmap(t(X))

### References
- [Swiss-Prot - Wikipedia](https://ja.wikipedia.org/wiki/Swiss-Prot)
- [Nucleic Acids Res. 2015 Jan;43(Database issue):D204-12. UniProt: a hub for protein information.](http://www.ncbi.nlm.nih.gov/pubmed/25348405)
- [Database (Oxford). 2014 Mar 12;2014:bau016. Expert curation in UniProtKB: a case study on dealing with conflicting and erroneous data.](http://www.ncbi.nlm.nih.gov/pubmed/24622611)


----------

## Silva rRNA database
[リボソームRNA](https://ja.wikipedia.org/wiki/リボソームRNA)データベース

### Website
[Silva rRNA database](https://www.arb-silva.de)の[Download → Archive](https://www.arb-silva.de/download/archive/)をクリックし、[current](https://www.arb-silva.de/no_cache/download/archive/current/)/[Exports](https://www.arb-silva.de/no_cache/download/archive/current/Exports/)を開く。
例えば、*README.txt*ファイルを右クリックし、「リンクのURLをコピー (Copy Link)」する。

	README.txt
	SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
	SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz.md5

### Download

[ターミナル](http://techacademy.jp/magazine/5155)を開き、`bash`を起動する:  

    bash

プロジェクト・ディレクトリを作成し移動する:  

    mkdir -p ~/projects/silva/data
    cd ~/projects/silva/data/

*README.txt*ファイル、塩基配列データの圧縮FASTAファイル（*.fasta.gz*）と[MD5](https://ja.wikipedia.org/wiki/MD5)[チェックサム](https://ja.wikipedia.org/wiki/チェックサム)ファイル（*.fasta.gz.md5*）を[`wget`](http://blog.layer8.sh/ja/2012/03/31/wget_command/)でダウンロードする:  

    wget --background \
     https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/README.txt \
     https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz \
     https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz.md5

または

    wget --background https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/{README.txt,SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz.md5,SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz}

`tail -f`でファイル出力を監視する（Control-Cで動作中のプロセスを停止）:  

    # Use `tail -f` to constantly monitor files (use Control-C to stop)
    tail -f wget-log

`md5`コマンドでチェックサムを計算し、公表されている値（492e513a8cc2de298a9b4121c9696278）と一致するか確認する:  

    md5 SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz
    cat SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz.md5

[Markdownノート（`README.md`）記載例](https://dl.dropboxusercontent.com/u/33495171/introBI/markdown/silva_rrna/README.md)

### Inspecting Data

`ls -l`でファイルの詳細情報を表示する:  

    ls -l

`gunzip`コマンドで解凍する:  

    gunzip -c SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta.gz > SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta

`ls -lh`でファイルサイズを確認する:  

    ls -lh

変数に値を割り当てる:  

    DB="SILVA_128_SSURef_Nr99_tax_silva_trunc.fasta"

`head`コマンドを用いて、ファイルの先頭部分を表示する:  

    head -n 3 $DB

`grep`コマンドを用いて、FASTA形式ファイルのヘッダ（`>`で始まる行）を抽出する:  

    grep "^>" $DB

[パイプ](https://ja.wikipedia.org/wiki/パイプ_%28コンピュータ%29)でプログラムの入出力をつなぐ。

配列の登録件数をカウントする:  

    grep "^>" $DB | wc -l

"Bacillus"にマッチする行を抽出する:  

    grep "^>" $DB | grep -c "Bacillus"

`grep`コマンドは、`-c`オプションでパターンにマッチした行数を表示し、`-i`オプションで大文字小文字を区別しない（ignore case）。

    grep "^>" $DB | grep -ci "Bacillus"

[FAQs](https://www.arb-silva.de/documentation/faqs/) "The taxonomic paths are standardized to six ranks; Domain, Phylum, Class, Order, Family and Genus."

Unixコマンド（`grep, cut, sort, uniq`）を組み合わせて、データの列を要約する。  
分類階級の[ドメイン](https://ja.wikipedia.org/wiki/ドメイン_%28分類学%29)：古細菌、真正細菌、真核生物（Domain: Archaea, Bacteria, Eukarya）をカウントする:  

    grep "^>" $DB | cut -d" " -f2 | cut -d";" -f1 | sort | uniq -c

ファイルで真正細菌（"Bacteria"）にマッチする行を抽出し、分類階級の[門](https://ja.wikipedia.org/wiki/門_%28分類学%29)（Phylum）をカウントし、`output.txt`ファイルへ出力する:  

    grep "^>" $DB | grep "Bacteria" | cut -d" " -f2 | cut -d";" -f2 | sort | uniq -c | sort -nr > output.txt

出力ファイルを確認する:  

    cat output.txt

### References
- [SILVA ribosomal RNA database - Wikipedia](https://en.wikipedia.org/wiki/SILVA_ribosomal_RNA_database)
- [blastn for silva database (SE)](http://cell-innovation.nig.ac.jp/wiki2/tiki-index.php?page=P000001306)

----------

## NCBI assembly summary
[NCBI](https://ja.wikipedia.org/wiki/国立生物工学情報センター)のゲノム配列のメタデータが記載されている。

### Website
[National Center for Biotechnology Information](http://www.ncbi.nlm.nih.gov)右下[NCBI FTP Site](ftp://ftp.ncbi.nlm.nih.gov/)を開き、[genomes](ftp://ftp.ncbi.nlm.nih.gov/genomes/)/[ASSEMBLY_REPORTS](ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/)に移動する。
 <ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt> で内容を確認:  

	assembly_summary_genbank.txt            - current GenBank genome assemblies
	assembly_summary_refseq.txt             - current RefSeq genome assemblies

### Download

[ターミナル](http://techacademy.jp/magazine/5155)を開き、`bash`を起動する:  

    bash

プロジェクト・ディレクトリを作成し移動する:  

    mkdir -p ~/projects/ncbi_assembly_summary/data
    cd ~/projects/ncbi_assembly_summary/data/

GenBankとRefSeqのゲノム配列のメタデータを記載したファイルを`wget`でダウンロードする:  

    # Download the two master assembly summary files that report assembly meta-data
    wget --background \
     ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt \
     ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

`tail -f`でファイル出力を監視する（Control-Cで動作中のプロセスを停止）:  

    # Use `tail -f` to constantly monitor files (use Control-C to stop)
    tail -f wget-log

ファイルのヘッダ（`#`で始まる行）を確認する:  

    grep "^#" assembly_summary_genbank.txt
    grep "^#" assembly_summary_refseq.txt

アセンブリの状況（assembly_level: Contig, Scaffold, Complete Genome）を確認する:  

    grep -v "^#" assembly_summary_genbank.txt | cut -f12 | sort | uniq -c
    grep -v "^#" assembly_summary_refseq.txt | cut -f12 | sort | uniq -c

[腸管出血性大腸菌O157](https://ja.wikipedia.org/wiki/O157) [Escherichia coli O157:H7 Sakai](http://integbio.jp/dbcatalog/record/nbdc00058) のゲノム配列のメタデータを確認し、RefSeq完全ゲノム('Complete Genome')配列データの最新版('latest')のURLを抽出する:  

    NAME="O157.*Sakai"
    grep "$NAME" assembly_summary_genbank.txt
    grep "$NAME" assembly_summary_refseq.txt
    ftpdirpaths=(`awk -F "\t" '$8 ~ /'"$NAME"'/ && $11=="latest" && $12 ~ /Complete Genome/ {print $20}' assembly_summary_refseq.txt`)

    # 抽出されたURLの数を確認する:  
    echo ${#ftpdirpaths[@]}

    # 抽出されたURLを表示する:  
    echo ${ftpdirpaths[@]}

抽出されたURL <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000008865.1_ASM886v1> をブラウザ（FirefoxまたはChrome）で開く。

	ファイル名	データファイルの内容
	*.gbff.gz	GenBank flat file format - GenBank形式ファイル
	*.fna.gz	FASTA Nucleic Acids - ゲノム塩基配列
	*.faa.gz	FASTA Amino Acids - 各タンパク質のアミノ酸配列

*README.txt*ファイル、[MD5](https://ja.wikipedia.org/wiki/MD5)[チェックサム](https://ja.wikipedia.org/wiki/チェックサム)ファイル（*md5checksums.txt*）、[GenaBank](http://bi.biopapyrus.net/biodb/genbank.html)形式と[FASTA](https://ja.wikipedia.org/wiki/FASTA)形式の圧縮ファイル（*.gz*）を`wget`でダウンロードする:  

    for URL in ${ftpdirpaths[@]}
    do
      wget $URL/{README.txt,md5checksums.txt,*.gbff.gz,*.fna.gz,*.faa.gz}    done

`md5`コマンドでチェックサムを計算し、公表されている値と一致するか確認する:  

    md5 *.gz
    cat md5checksums.txt

[Markdownノート（`README.md`）記載例](https://dl.dropboxusercontent.com/u/33495171/introBI/markdown/ncbi_assembly_summary/README.md)

### References
- [RefSeq - JI 井上 潤](http://www.geocities.jp/ancientfishtree/RefSeq.html)
- [RefSeq | 重複のない生物の遺伝子データベース（ゲノムデータベース）](http://bi.biopapyrus.net/biodb/refseq.html)
- [What is the difference between RefSeq and GenBank?](https://www.ncbi.nlm.nih.gov/books/NBK50679/#RefSeqFAQ.what_is_the_difference_between_1)
- [Genomes Download FAQ](http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/)
 - [1. What is the best protocol to use to download data?](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#protocols)
 - [11. How can I find the sequence and annotation of my genome of interest?](http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#howtofind)
 - [14. How can I download RefSeq data for all complete bacterial genomes?](http://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete)

----------







