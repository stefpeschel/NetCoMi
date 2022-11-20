
## NetCoMi 1.1.0 <img src="man/figures/NetCoMi_logo_800x400_300dpi.png" align="right" width="200" />

### New features

- **renameTaxa()**: New function for renaming taxa in a taxonomic table.
  It comes with functionality for making unknown and unclassified taxa
  unique and substituting them by the next higher known taxonomic level.
  E.g., an unknown genus “g\_\_“, where family is the next higher known
  level, can automatically be renamed to”1_Streptococcaceae(F)“.
  User-defined patterns determine the format of known and substituted
  names. Unknown names (e.g., NAs) and unclassified taxa can be handled
  separately. Duplicated names within one or more chosen ranks can also
  be made unique by numbering them consecutively.

- **editLabels()**: New function for editing node labels, i.e.,
  shortening to a certain length and removing unwanted characters. It is
  used by NetCoMi’s plot functions plot.microNetProps() and
  plot.diffnet().

- In `netCompare()`: The **adjusted Rand index** is also computed for
  the largest connected component (LCC). The summary method has been
  adapted.

- Argument **“testRand”** added to `netCompare()`. Performing a
  permutation test for the adjusted Rand index can now be disabled to
  save run time.

- **Graphlet-based network measures** implemented. NetCoMi contains two
  new exported functions **`calcGCM()`** and **`calcGCD()`** to compute
  the Graphlet Correlation Matrix (GCM) of a network and the Graphlet
  Correlation Distance (GCD) between two networks. **Orbits** for
  graphlets with up to four nodes are considered. Furthermore, the GCM
  is computed with `netAnalyze()` and the GCD with `netCompare()` (for
  the whole network and the largest connected component, respectively).
  Also the orbit counts are returned. The GCD is added to the summary
  for class `microNetComp` objects returned by `netCompare()`.

- **Significance test for the GCD**: If permutation tests are conducted
  with `netCompare()`, the GCD is tested for being significantly
  different from zero.

- New function **`testGCM()`** to **test graphlet-based measures** for
  significance. For a single GCM, the correlations are tested for being
  significantly different from zero. If two GCMs are given, it is tested
  if the correlations are significantly different between the two
  groups, that is, the absolute differences between correlations (
  $|gc1_{ij}-gc2_{ij}|$ ) are tested for being different from zero.

- New function **`plotHeat()`** for plotting a mixed heatmap where, for
  instance, values are shown in the upper triangle and corresponding
  p-values or significance codes in the lower triangle. The function is
  used for plotting heatmaps of the GCMs, but could also be used for
  association matrices.

- `netAnalyze()` now by default returns a **heatmap of the GCM(s)** with
  graphlet correlations in the upper triangle and significance codes in
  the lower triangle.

- Argument **“doPlot”** added to `plot.microNetProps()` to suppress the
  plot if only the return value is of interest.

- New **“show”** arguments are added to the summary methods for class
  `microNetProps` and `microNetComp` objects. They specify which network
  properties should be printed in the summary. See the help pages of
  `summary.microNetProps` and `summary.microNetComp()` for details.

- New **zero replacement** method **“pseudoZO”** available in
  `netConstruct()`. Instead of adding the desired pseudo count to the
  whole count matrix, it is added to zero counts only if `pseudoZO` is
  chosen. The behavior of “pseudo” (a further available method where a
  pseudo count is added to all counts) has not changed. Adding a pseudo
  count only to zeros preserves the ratios between non-zero counts,
  which is desirable.

- `createAssoPerm()` now accepts objects of class `microNet` as input
  (in addition to objects of class `microNetProps`).

- **`SPRING's`** fast version of latent correlation computation
  (implemented in [mixedCCA](https://github.com/irinagain/mixedCCA)) is
  available again. It can be used by setting the `netConstruct()`
  parameter `measurePar$Rmethod` to “approx”, which is now the default
  again.

- The function **`multAdjust()`** now has an argument `pTrueNull` to
  pre-define the proportion of true null hypotheses for the adaptive BH
  method.

- `netConstruct()` has a new argument **`assoBoot`**, which enables the
  computation of bootstrap association matrices outside netConstruct()
  if **bootstrapping** is used for sparsification. An example has been
  added to the help page `?netConstruct`. This feature might be useful
  for very large association matrices (for which the working memory
  might reach its limit).

### Bug fixes

- In `netConstruct()`:

  - Using **“bootstrap”** as sparsification method in combination with
    one of the association methods “bicor”, “cclasso”, “ccrepe”, or
    “gcoda” led to the error:
    `argument "verbose" is missing, with no default`, which has been
    fixed.
  - The **“signedPos”** transformation did not work properly.
    Dissimilarities corresponding to negative correlations were set to
    zero instead of infinity.

- In `editLabels()`: The function (and thus also `plot.microNetProps`)
  threw an error if taxa have been renamed with `renameTaxa` and the
  data contain more than 9 taxa with equal names, so that double-digit
  numbers were added to avoid duplicates.

- Issues in network analysis and plotting if association matrices are
  used for network construction, but **row and/or column names are
  missing**. (issue
  [\#65](https://github.com/stefpeschel/NetCoMi/issues/65))

- `diffnet()` threw an error if association matrices are used for
  network construction instead of count matrices. (issue
  [\#66](https://github.com/stefpeschel/NetCoMi/issues/66))

- In `plot.microNetProps()`:

  - The function now directly returns an error if `x` has not the
    expected class.
  - The `cut` parameter could not be changed.

- In **`cclasso()`**: In rare cases, the function produced complex
  numbers, which led to an error.

### Further changes

- In **permutation tests**: The permuted group labels must now be
  different from the original group vector. In other words, the original
  group vector is strictly avoided in the matrix with permuted group
  labels. So far, only duplicates were avoided. Only in exact
  permutation tests (if `nPerm` equals the possible number of
  permutations), the original group vector is still included in the
  permutation matrix. The calculation of p-values has been adapted to
  the new behavior: *p=B/N* for exact p-values and *p=(B+1)/(N+1)* for
  approximated p-values, where *B* is the number of permutation test
  statistics being larger than or equal to the observed one, and *N* is
  the number of permutations. So far, *p=(B+1)/(N+1)* has been used in
  all cases.

- In `plot.microNetProps()`:

  - The default of `shortenLabels` is now “none”, i.e. the **labels are
    not shortened by default**, to avoid confusion about the node
    labels.
  - The **edge filter** (specified via `edgeFilter` and
    `edgeInvisFilter`) now refers to the estimated
    association/dissimilarities instead of edge weights. E.g., setting
    the threshold to 0.3 for an association network hides edges with a
    corresponding absolute association below 0.3 even though the edge
    weight might be different (depending on the transformation used for
    network construction). (issue
    [\#26](https://github.com/stefpeschel/NetCoMi/issues/26))
  - If two networks are constructed and the **`cut`** parameter is not
    user-defined, the mean of the two determined cut parameters is now
    used for both networks so that edge thicknesses are comparable.

- More expressive messages and errors in `diffnet` and `plot.diffnet` if
  no **differential associations** are detected.

- New function **`.suppress_warnings()`** to suppress certain warnings
  returned by external functions.

- In `netConstruct` if **“multRepl”** is used for zero handling: The
  warning about the proportion of zeros is suppressed by setting the
  `multRepl()` parameter “z.warning” to 1.

- The functions **`makeCluster`** and **`stopCluster`** from `parallel`
  package are now used for parallel computation because those from
  `snow` package sometimes led to problems on Unix machines.

### Style

- The whole R code has been reformatted to follow general conventions.

- The element `"clustering_lcc"` as part of the `netAnalyze` output has
  changed to `"clusteringLCC"` to be in line with the remaining output.

- Input argument checking of exported function has been revised. New
  functions `.checkArgsXxx()` are added to perform argument checking
  outside the main functions.

- Non-exported functions have been renamed to follow general naming
  conventions, i.e. that of
  [Bioconductor](https://contributions.bioconductor.org/r-code.html):

  - Use camelCase for all functions.
  - Non-exported functions have prefix “.”
  - The following functions have been renamed:

| Old names              | New names             |
|:-----------------------|:----------------------|
| boottest               | .boottest             |
| calc_association       | .calcAssociation      |
| calc_diff_props        | .calcDiffProps        |
| calc_jaccard           | .calcJaccard          |
| calc_props             | .calcProps            |
| diff_connect_pairs     | .diffConnectPairs     |
| diff_connect_variables | .diffConnectVariables |
| diff_connect_network   | .diffConnectNetwork   |
| filter_edges           | .filterEdges          |
| filter_nodes           | .filterNodes          |
| filter_samples         | .filterSamples        |
| filter_taxa            | .filterTaxa           |
| first_unequal_element  | .firstUnequalElement  |
| get_clust_cols         | .getClustCols         |
| get_node_size          | .getNodeSize          |
| get_perm_group_mat     | .getPermGroupMat      |
| get_vec_names          | .getVecNames          |
| norm_counts            | .normCounts           |
| permtest_diff_asso     | .permTestDiffAsso     |
| scale_diss             | .scaleDiss            |
| sparsify               | .sparsify             |
| trans_to_diss          | .transToDiss          |
| trans_to_sim           | .transToSim           |
| trans_to_adja          | .transToAdja          |
| zero_treat             | .zeroTreat            |

## NetCoMi 1.0.3

This is a minor release with some bug fixes and changes in the
documentation.

### Bug fixes

- `netConstruct()` threw an error if the data had no row and/or column
  names, which is fixed.

- An edge list is added to the output of `netConstruct()` (issue
  [\#41](https://github.com/stefpeschel/NetCoMi/issues/41)). See the
  help page for details.

- `SPRING`’s fast version of latent correlation computation (implemented
  in [mixedCCA](https://github.com/irinagain/mixedCCA)) is currently not
  available due to deprecation of the R package `chebpol`. The issue is
  fixed by setting the `netConstruct()` parameter `measurePar$Rmethod`
  internally to “original” if SPRING is used for association estimation.

- In `plot.microNetProps()`: The `xpd` parameter is changed to `NA` so
  that plotting outside the plot region is possible (useful for legends
  or additional text).

- Labels in the network plot can now be suppressed by setting
  `labels = FALSE` (issue
  [\#43](https://github.com/stefpeschel/NetCoMi/issues/43))

- The `netCompare()` function threw an error if one of the permutation
  networks was empty, i.e. had no edges with weight different from zero
  (issue [\#38](https://github.com/stefpeschel/NetCoMi/issues/38)),
  which is now fixed.

- Fix issues [\#29](https://github.com/stefpeschel/NetCoMi/issues/29)
  and [\#40](https://github.com/stefpeschel/NetCoMi/issues/40), where
  permutation tests did not terminate for small sample sizes. Now, if
  the possible number of permutations (resulting from the sample size)
  is smaller than that defined by the user, the function stops and
  returns an error.

- Fix a bug in `diffnet()` (issue
  [\#51](https://github.com/stefpeschel/NetCoMi/issues/51)), where
  colors in differential networks could not be changed.

- `diffnet()` threw an error if the `netConstruct()` argument
  `jointPrepro` was set to `TRUE`.

## NetCoMi 1.0.2

This release includes a range of new features and fixes known bugs and
issues.

### New features

#### Improved installation process

Packages that are optionally required in certain settings are not
installed together with `NetCoMi` anymore. Instead, there is a new
function `installNetCoMiPacks()` for installing the remaining packages.
If not installed via `installNetCoMiPacks()`, the required package is
installed by the respective NetCoMi function when needed.

#### installNetCoMiPacks()

New function for installing the R packages used in NetCoMi not listed as
`dependencies` or `imports` in NetCoMi’s description file.

#### netConstruct()

- New argument `matchDesign`: Implements matched-group
  (i.e. matched-pair) designs, which are used for permutation tests in
  `netCompare()` and `diffnet()`. `c(1,2)`, for instance, means that one
  sample in the first group is matched to two samples in the second
  group. If the argument is not `NULL`, the matched-group design is kept
  when generating permuted data.

- New argument `jointPrepro`: Specifies whether two data sets (of group
  one and two) should be preprocessed together. Preprocessing includes
  sample and taxa filtering, zero treatment, and normalization. Defaults
  to `TRUE` if `data` and `group` are given, and to `FALSE` if `data`
  and `data2` are given, which is similar to the behavior of
  `NetCoMi 1.0.1`. For dissimilarity networks, no joint preprocessing is
  possible.

- `mclr(){SPRING}` is now available as normalization method.

- `clr{SpiecEasi}` is used for centered log-ratio transformation instead
  of `cenLR(){robCompositions}`.

- `"symBetaMode"` is accepted as list element of `measurePar`, which is
  passed to `symBeta(){SpiecEasi}`. Only needed for SpiecEasi or SPRING
  associations.

- The pseudocount (if `zeroMethod = "pseudo"`) may be freely specified.
  In v1.0.1, only unit pseudocounts were possible.

#### netAnalyze()

- Global network properties are now computed for the whole network as
  well as for the largest connected component (LCC). The summary of
  network properties now contains for the whole network only statistics
  that are not based on shortest paths (or, more generally, also
  meaningful for disconnected networks). For the LCC, all global
  properties available in NetCoMi are shown.

- New global network properties (see the docu of `netAnalyze()` for
  definitions):

  - Number of components (only whole network)
  - Relative LCC size (only LCC)
  - Positive edge percentage
  - Natural connectivity
  - Average dissimilarity (only meaningful for the LCC)
  - Average path length (only meaningful for the LCC)

- New argument `centrLCC`: Specifies whether to compute centralities
  only for the LCC. If `TRUE`, centrality values of disconnected
  components are zero.

- New argument `avDissIgnoreInf`: Indicates whether infinite values
  should be ignored in the average dissimilarity. If `FALSE`, infinities
  are set to 1.

- New argument `sPathAlgo`: Algorithm used for computing shortest paths

- New argument `sPathNorm`: Indicates whether shortest paths should be
  normalized by average dissimilarity to improve interpretability.

- New argument `normNatConnect`: Indicates whether to normalize natural
  connectivity values.

- New argument `weightClustCoef`: Specifies the algorithm used for
  computing the global clustering coefficient. If `FALSE`,
  `transitivity(){igraph}` with `type = "global"` is used (similar to
  `NetCoMi 1.0.1`). If `TRUE`, the local clustering coefficient is
  computed using `transitivity(){igraph}` with `type = "barrat"`. The
  global clustering coefficient is then the arithmetic mean of local
  values.

- Argument `connect` has been changed to `connectivity`.

- Documentation extended by definitions of network properties.

#### summary.microNetProps()

- New argument `clusterLCC`: Indicates whether clusters should be shown
  for the whole network or only for the LCC.

- The `print` method for `summary.microNetProps` was completely revised.

#### plot.microNetProps()

- All normalization methods available for network construction can now
  be used for scaling node sizes (argument `nodeSize`).

- New argument `normPar`: Optional parameters used for normalization.

- Usage of `colorVec` changed: Node colors can now be set separately in
  both groups (`colorVec` can be a single vector or a list with two
  vectors). Usage depends on `nodeColor` (see docu of `colorVec`).

- New argument `sameFeatCol`: If `nodeColor = "feature"` and `colorVec`
  is not given, `sameFeatCol` indicates whether same features should
  have same colors in both groups.

- Argument `colorNegAsso` has been renamed to `negDiffCol`. Using the
  old name leads to a warning.

- New functionality for using the same layout in both groups (if two
  networks are plotted). In addition to computing the layout for one
  group and adopting it for the other group, a union of both layout can
  be computed and used in both groups so that nodes are placed as
  optimal as possible equally for both networks. This option is applied
  via `sameLayout = TRUE` and `layoutGroup = "union"`. Many thanks to
  [Christian L. Müller](https://github.com/muellsen?tab=followers) and
  [Alice Sommer](https://www.iq.harvard.edu/people/alice-sommer) for
  providing the idea and R code for this new feature!

#### netCompare()

- New arguments for storing association and count matrices of the
  permuted data into an external file:

  - `fileLoadAssoPerm`
  - `fileLoadCountsPerm`
  - `storeAssoPerm`
  - `fileStoreAssoPerm`
  - `storeCountsPerm`
  - `fileStoreCountsPerm`

- New argument `returnPermProps`: If `TRUE`, global network properties
  and the respective absolute group differences of the permuted data are
  returned.

- New argument `returnPermCentr`: If `TRUE`, the computed centrality
  values and the respective absolute group differences of the permuted
  data are returned as list with a matrix for each centrality measure.

- The arguments `assoPerm` and `dissPerm` are still existent for
  compatibility with `NetCoMi 1.0.1` but the former elements `assoPerm`
  and `dissPerm` are not returned anymore (matrices are stored in an
  external file instead).

#### createAssoPerm()

New function for creating association/dissimilarity matrices for
permuted count data. The stored count or association/dissimilarity
matrices can then be passed to `netCompare()` or `diffnet()` to decrease
runtime. The function also allows to generate a matrix permuted group
labels without computing associations. Using this matrix,
`createAssoPerm()` furthermore allows to estimate the permutation
associations/dissimilarities in blocks (by passing only a subset of the
permuted group matrix to `createAssoPerm()`).

#### summary.microNetComp()

Summary method has been adapted to the new network properties (analogous
to the summary of `microNetProps` objects, which are returned from
`netAnalyze()`)

#### diffnet()

- New arguments for storing association and count matrices of the
  permuted data into an external file:

  - `fileLoadAssoPerm`
  - `fileLoadCountsPerm`
  - `storeAssoPerm`
  - `fileStoreAssoPerm`
  - `storeCountsPerm`
  - `fileStoreCountsPerm`

- The argument `assoPerm` is still existent for compatibility with
  `NetCoMi 1.0.1` but the former element `assoPerm` is not returned
  anymore (matrices are stored in an external file instead).

- Changed output: For permutation tests and Fisher’s z-test, a vector
  and matrix with p-values and the corresponding matrix with group
  differences are returned for both with and without multiple testing
  adjustment.

- Documentation has been revised.

#### plot.diffnet()

- New argument `adjusted`: Indicates whether the adjacency matrix
  (matrix with group differences) based on adjusted or unadjusted
  p-values should be plotted.

- New argument `legendPos` for positioning the legend.

- New argument `legendArgs` for specifying further arguments passed to
  `legend`.

#### colToTransp()

- The function is now exported and its name has changed from
  `col_to_transp()` to `colToTransp()`. The function expects a color
  vector as input and adds transparency to each color.

### Bug fixes

The major issues fixed in this release are:

- The following error is solved:
  `Error in update.list(...): argument "new" is  missing`. The error was
  caused by a conflict between `SpiecEasi` and `metagenomeSeq`, in
  particular by `gplot` as a dependency of `metagenomeSeq`. A former
  version of `gplot` was dependend on `gdata`, which caused the
  conflict. So, please update `gplot` and remove the package `gdata` to
  fix the error.

- `sparcc()` from SpiecEasi package is now used for estimating SparCC
  associations. For some users, NetCoMi’s `Rccp` implementation of
  SparCC caused errors when installing NetCoMi. If these are fixed, the
  Rcpp implementation will be included again, so that users can decide
  between the two SparCC versions.

- VST transformations are now computed correctly.

- Error when plotting two networks, where one network is empty, has been
  fixed.
