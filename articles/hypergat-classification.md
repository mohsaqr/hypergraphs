# Classifying documents with existing subject labels

A document classifier learns from documents whose subjects have already
been assigned. Its usefulness is measured by the subjects it predicts
for other documents with known labels. The R8 news corpus has eight
existing subject labels, including earnings, acquisitions and trade, and
a fixed training and test split. HyperGAT learns from the training
articles and classifies the test articles with the same fitted network.

## Existing labels

R8 is a subset of Reuters newswire articles with one of eight subject
labels. The [corpus
text](https://raw.githubusercontent.com/kaize0409/HyperGAT_TextClassification/9af29b6837b95344bbc5a567e0673b41c616c88f/data/R8_corpus.txt)
and [label
file](https://raw.githubusercontent.com/kaize0409/HyperGAT_TextClassification/9af29b6837b95344bbc5a567e0673b41c616c88f/data/R8_labels.txt)
come from a fixed revision of the official HyperGAT implementation.
Download both files into the working directory before running the code.
Fitting requires the `torch` R package and its runtime.

The tables and figures below retain a fitted reference run with
`seed = 1`. Package builds use these saved results. The displayed code
reproduces the analysis; numerical results can vary across torch
versions and hardware.

The files retain sentence punctuation and the original subject labels.
The training and test articles are kept separate. The test labels are
used to evaluate predictions after fitting.

``` r

text <- enc2utf8(readLines(con = "R8_corpus.txt", encoding = "latin1",
                          warn = FALSE))
meta <- read.delim(file = "R8_labels.txt", header = FALSE,
                   col.names = c("row", "split", "label"),
                   colClasses = "character", quote = "")
articles <- transform(`_data` = meta,
                       node = sprintf("R8_%05d", seq_along(text)), text = text)
training <- subset(x = articles, subset = split == "train",
                   select = c(node, text, label))
testing <- subset(x = articles, subset = split == "test",
                  select = c(node, text, label))
```

``` r

xtabs(formula = ~ label + split, data = articles)
```

              split
    label      test train
      acq       696  1596
      crude     121   253
      earn     1083  2840
      grain      10    41
      interest   81   191
      money-fx   87   222
      ship       36   108
      trade      76   250

There are 5501 training articles and 2190 test articles. Earnings
accounts for 51.6% of the training set. A classifier that always chooses
the most common training subject provides a simple comparison for
overall accuracy. Balanced accuracy also measures how well the smaller
subjects are recovered.

## Fit once

HyperGAT represents each article as a hypergraph. Its nodes are the
unique words, and each sentence is a hyperedge containing its words. The
network learns word representations and aggregates them through these
hyperedges to predict a subject. The default construction uses
sentences. Topic hyperedges are optional.

``` r

fit <- hg_hypergat(x = training, column = "text", id = "node",
                   labels = "label", min_count = 5, seed = 1)
print(x = fit)
```

    HyperGAT classifier: 5501 documents, 8 classes, 5501 labelled
    No held-out evaluation. predict() reuses the fitted model.

The fitted object retains its network and vocabulary. The vocabulary
comes from the known training labels. The default validation split
selects the training epoch; it supplies no estimate of performance on
the test set. For a corpus without a supplied test split,
`holdout = 0.2` provides an internal evaluation using a stratified share
of the known labels.

``` r

plot(x = fit, type = "history")
```

![Training loss of the fitted document
classifier](../../../../../../tmp/RtmpBBTYOU/temp_libpath18ba3e3e96de/hypergraphs/extdata/hypergat-classification/history.png)

## Predict independent documents

[`predict()`](https://rdrr.io/r/stats/predict.html) applies the fitted
network to the test articles. Supplying their existing label column
evaluates the predictions. These labels do not update the network or
vocabulary.

``` r

evaluation <- predict(object = fit, newdata = testing, labels = "label")
print(x = evaluation)
```

    Classification (hg_hypergat prediction): 2190 documents, 2190 held out (100% of the labelled)
    Held-out accuracy 0.963, balanced accuracy 0.898 (chance 0.125)

        class n_test correct    recall precision
          acq    696     678 0.9741379 0.9854651
        crude    121     113 0.9338843 0.8968254
         earn   1083    1066 0.9843029 0.9897864
        grain     10       9 0.9000000 0.9000000
     interest     81      65 0.8024691 0.8552632
     money-fx     87      78 0.8965517 0.8297872
         ship     36      26 0.7222222 0.7428571
        trade     76      74 0.9736842 0.9024390

The accuracy is 96.3% on all 2190 test articles, compared with 49.5% for
always choosing earn. Balanced accuracy is 89.8%. The model can score
2188 of the test articles. Articles with no word in the fitted
vocabulary retain their rows and count as errors. These results measure
recovery of the eight existing subjects in this test set.

``` r

plot(x = evaluation)
```

![Existing subject against predicted subject on the independent test
set](../../../../../../tmp/RtmpBBTYOU/temp_libpath18ba3e3e96de/hypergraphs/extdata/hypergat-classification/confusion.png)

## Read the errors

The confusion plot shows which subjects are confused. Precision and
recall give the performance for each subject. The document table returns
the text together with its observed and predicted subject, so an error
can be inspected directly.

``` r

hg_get(x = evaluation, what = "classes")
```

| class    | n_test | correct |    recall | precision |
|:---------|-------:|--------:|----------:|----------:|
| acq      |    696 |     678 | 0.9741379 | 0.9854651 |
| crude    |    121 |     113 | 0.9338843 | 0.8968254 |
| earn     |   1083 |    1066 | 0.9843029 | 0.9897864 |
| grain    |     10 |       9 | 0.9000000 | 0.9000000 |
| interest |     81 |      65 | 0.8024691 | 0.8552632 |
| money-fx |     87 |      78 | 0.8965517 | 0.8297872 |
| ship     |     36 |      26 | 0.7222222 | 0.7428571 |
| trade    |     76 |      74 | 0.9736842 | 0.9024390 |

``` r

hg_get(x = evaluation, what = "documents", correct = FALSE, top = 1)
```

| node | label | predicted | score | margin | split | correct | text |
|:---|:---|:---|---:|---:|:---|:---|:---|
| R8_05513 | money-fx | interest | 0.8790182 | 0.7759423 | test | FALSE | LAWSON SEES NO CHANGE IN U.K. MONETARY POLICY. British |

Chancellor of the Exchequer Nigel Lawson said he saw no immediate
implications for British monetary policy arising from the Group of Seven
meeting yesterday. “Exchange rate stability is in the U.K.’s interest,”
he told journalists. Asked what it meant for U.K. monetary policy, he
said, “No, I do not think there are any immediate implications.” \|
\|R8_05601 \|trade \|crude \| 0.4464583\| 0.0061157\|test \|FALSE
\|SPRINKEL SAYS TAX HIKE WOULD NOT REDUCE DEFICIT. Council of Economic
Advisers chairman Beryl Sprinkel said the Reagan Administration remains
strongly opposed to a tax increase, including 18 billion dlrs of new
revenues in the budget plan by Congressional Democrats. “We believe that
significant increases in taxes would not reduce deficits and could have
adverse effects on growth,” Sprinkel told the House Rules Committee. He
said the Administration wanted to continue its policy of gradually
reducing deficits through restraining government spending and promoting
economic growth. Sprinkel said cutting the budget deficit was the best
way to lower the trade deficit. \| \|R8_05668 \|ship \|earn \|
0.4305401\| 0.2379675\|test \|FALSE \|TODD SHIPYARDS STRUCK ON WEST
COAST. Todd Shipyards Corp said production workers represented by the
multi-union Pacific Coast Metal Trades District Council at its San
Francisco division struck on April Six. It said negotiations are
expected to resume at the end of this month. Todd also said the
collective bargaining division in effect at its Galveston Division
expires April 17, and negotiations with the Galveston Metal Trades
Council are continuing. The company said results of balloting on a new
collective bargaining agreement proposal in its Seattle Division are
expected to be tabulated at the close of business tomorrow. The Pacific
Coast Council has recommended acceptance of that proposal by membership,
Todd said. \| \|R8_05701 \|interest \|money-fx \| 0.6417318\|
0.2867234\|test \|FALSE \|FED EXPECTED TO SET CUSTOMER REPURCHASES. The
Federal Reserve is expected to intervene in the government securities
market to supply temporary reserves indirectly via customer repurchase
agreements, economists said. Economists expect the Fed to execute
2.0-2.5 billion dlrs of customer repos to offset pressures from the end
of the two-week bank reserve maintenance period today. Some also look
for a permanent reserve injection to offset seasonal pressures via an
outright purchase of bills or coupons this afternoon. The Federal funds
rate opened at 6-3/8 pct and remained at that level, up from yesterday’s
6.17 pct average. \| \|R8_05795 \|crude \|acq \| 0.5709212\|
0.1564069\|test \|FALSE \|GULF CANADA ACQUIRES SUEZ OIL STAKE. Gulf
Canada Corp said it acquired a 25 pct working interest in the Gulf of
Suez oil concession for undisclosed terms. The company said its
agreement with operator Conoco Hurghada Inc and Hispanoil covered the
168,374-acre East Hurghada offshore concession. It said a 15.6 mln U.S.
dlr four-well program was planned for 1987. After the acquisition, which
is subject to Egyptian government approval, working interests in the
Hurghada block will be Conoco Hurghada at 45 pct, Hispanoil 30 pct and
Gulf Canada the balance. \| \|R8_05930 \|earn \|interest \| 0.7923861\|
0.6172475\|test \|FALSE \|BANKERS TRUST PUTS BRAZIL ON NON-ACCRUAL.
Bankers Trust New York Corp said it has placed its approximately 540 mln
dlrs of medium- and long-term loans to Brazil on non-accrual status and
that first-quarter net income will be reduced by about seven mln dlrs as
a result. Brazil suspended interest payments on its 68 billion dlrs of
medium- and long-term debt on February 22. U.S. banking regulations do
not require banks to stop accruing interest on loans until payments are
90 days overdue, but Bankers Trust said it acted now because of “the
high potential of a continued suspension that would result in reaching
the 90-day limit in the second quarter of 1987.” Assuming no cash
payments at current interest rates are received for the rest of 1987,
Bankers Trust estimated that full-year net income would be reduced by
about 30 mln dlrs. Bankers Trust said it assumes that debt negotiations
between Brazil and its commercial bank lenders will lead to the
resumption of interest payments. The negotiations resume in New York on
Friday when central bank governor Francisco Gros is expected to ask
banks for a 90-day rollover of some 9.5 billion dlrs of term debt that
matures on April 15. \| \|R8_07140 \|acq \|earn \| 0.6432334\|
0.3158012\|test \|FALSE \|CANADA STOCKS/DOME PETROLEUM LTD . Dome
Petroleum Ltd shares moved higher in the U.S. and Canada after
TransCanada PipeLines Ltd made a 4.3 billion Canadian dlr bid for Dome
and Dome said it is in talks with two other unidentified companies.
Market speculation is that the other two potential bidders are not
Canadian companies and DuPont’s

Conoco and Atlantic Richfield Co are mentioned as possibilities, Wilf
Gobert of Peters and Co Ltd said. Dome rose 1/4 to 1-1/8 on the American
Stock Exchange. TransCanada PipeLines was down 1/4 at 15-3/4 on the New
York Stock Exchange. Dome was the most active stock on the Toronto
exchange at 1.50 dlrs per share, up 37 cts. Gobert characterized the
market action in Dome as “awfully optimistic” but said investors are
hoping for a competing offer to the shareholders. TransCanada PipeLines’
offer is to Dome management, not to shareholders. However, it proposes
issuing new equity in a subsidiary that would operate Dome assets.
Current Dome shareholders would own 20 pct of the new subsidiary. \|
\|R8_07691 \|grain \|interest \| 0.9186397\| 0.8909166\|test \|FALSE
\|U.S. SENATE PANEL VOTES TO LIMIT COUNTY LOAN RATE CHANGES STARTING
WITH 1988 CROPS. U.S. SENATE PANEL VOTES TO LIMIT COUNTY LOAN RATE
CHANGES STARTING WITH 1988 CROPS \|

A new article without a label is classified with
`predict(object = fit, newdata = new_articles)`. The output contains its
predicted subject, winning softmax score and margin over the next
subject. The score is a model output rather than a calibrated
probability that the subject is correct. Every scorable article receives
one of the fitted subjects. These labels are appropriate for documents
belonging to this subject system; this evaluation does not establish
accuracy for other corpora or subjects.

## Attention as a diagnostic

Attention weights determine internal aggregation. Word attention sums to
one over the words of a hyperedge. Edge attention sums to one over the
hyperedges containing a word. A word occurring in only one hyperedge
therefore assigns that hyperedge a weight of one automatically. The mean
of these edge weights can be high because many words are exclusive to a
sentence.

The diagnostic table includes the mean uniform weight implied by the
same memberships and the fraction of exclusive words. These values make
the normalization visible. The illustrated article is a short earnings
report with a headline and reported results. Its sentences have distinct
word sets, so the two hyperedges and their overlap can be read. It is
selected for readability, independently of confidence or correctness.

``` r

first_article <- "R8_01282"
```

The hypergraph plot shows the words and the sentence hyperedges that the
network consumes for this article.

``` r

plot(x = fit, type = "hyperedges", node = first_article, edge_labels = TRUE,
     label_size = 3, titles = FALSE)
```

![Words joined by the sentence hyperedges of one training
article](../../../../../../tmp/RtmpBBTYOU/temp_libpath18ba3e3e96de/hypergraphs/extdata/hypergat-classification/hyperedges.png)

``` r

hg_get(x = fit, what = "hyperedges", node = first_article)
```

| node | label | predicted | hyperedge | kind | text | weight | top_word | n_words | uniform_weight | exclusive_fraction |
|:---|:---|:---|:---|:---|:---|---:|:---|---:|---:|---:|
| R8_01282 | earn | earn | sentence 2 | sentence | Shr 2-1/5 cts vs nil Net 156,726 vs 11,989 Sales |  |  |  |  |  |
| 1,157,883 vs 890,138 | 0.9166800 | nil | 6 | 0.9166667 | 0.8333333 |  |  |  |  |  |
| R8_01282 | earn | earn | sentence 1 | sentence | CCR VIDEO 1ST QTR NOV 30 NET | 0.9166533 | qtr | 6 | 0.9166667 | 0.8333333 |

``` r

plot(x = fit, type = "attention", node = first_article)
```

![Internal edge attention with its uniform normalization baseline for
one training
article](../../../../../../tmp/RtmpBBTYOU/temp_libpath18ba3e3e96de/hypergraphs/extdata/hypergat-classification/attention.png)

Crosses show the uniform normalization baseline and bars show the mean
learned edge attention. These internal weights do not measure a
sentence’s contribution to the predicted subject. Classification errors
are reviewed from the document text and the evaluation tables.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification.
*Proceedings of EMNLP 2020*, 4927–4936.
[Paper](https://aclanthology.org/2020.emnlp-main.399/). The
sentence-preserving corpus and existing labels are distributed in the
[official
implementation](https://github.com/kaize0409/HyperGAT_TextClassification).
