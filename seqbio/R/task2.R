#' @title Muestra la significaci√≥n estad√≠stica de un alineamiento por parejas.
#' @description Dados dos secuencias, hace shuffling en una de ellas y dibujamos
#' un histograma con el numero de scores, adem√°s de devolver una descripci√≥n.
#' @param seq1 Una secuencia en formato fasta
#' @param seq2 Otra secuencia en formato fasta
#' @param tipos Tipo de secuencia(Proteina o ADN)
#' @param tipoal Tipo de alineamiento(‚Äúlocal‚Äù o ‚Äúglobal‚Äù)
#' @param M Matriz de substituci√≥n(PAMn, BLOSUMn,...)
#' @param puntugap Gap: Open penalty, extended penalty
#' @param replcas N√∫mero de replicas
#' @param shuf Secuencia donde se hace shuffling: 1 o 2
#'
#' @return el histograma con los scores
#'
#' @examples
#'
#' seq_alignment(seq1 = "P0DP27.fa", seq2 = "Q9N1R0.fa", seq_type ="protein",
#'               seq_align = "local", mat = data("BLOSUM62"), gap = c(-12, -2), N = 100, shuff = 1)
#'
#' mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
#' seq_alignment(seq1 = "gi32141095_N_0.fa", seq2 = "gi32141095_N_1.fa", seq_type ="dna",
#'               seq_align = "global", mat = mat, gap = c(-8, -1), N = 2000, shuff)
#'
#' @references \url{https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter4.html}
#'
#'
#' @export
#


comparaseq <- function(seq1,seq2,tipos,tipoal,M, openpenalty, extendedpenalty,N,shuf){
  if(!require('Biostrings',quietly = T)){BiocManager::install("Biostrings")}
  require('Biostrings', quietly = T)

  #  ComprobaciÛn de los inputs:
  stopifnot(class(seq1) == "character", class(seq2) == "character",
             class(gap) == "numeric", class(N) == "numeric", class(shuff) == "numeric",
            length(gap) == 2, length(N) == 1)

  #  Lectura de las secuencias:

  if(seq_type == "dna"){
    s1 = readDNAStringSet(seq1, "fasta")
    s2 = readDNAStringSet(seq2, "fasta")
    name_mat = "DNAScoringMat"
  } else{
    s1 = readAAStringSet(seq1, "fasta")
    s2 = readAAStringSet(seq2, "fasta")
    name_mat = mat
  }

  #  Alineamiento original (”ptimo):
  S <- pairwiseAlignment(s1, s2, type = seq_align, substitutionMatrix = mat,
                         gapOpening = gap[1], gapExtension = gap[2], scoreOnly = TRUE)
  rep <- rep(0,N)
  if (shuf==1){
    seqaleat <- shuffling(as.character(s1[[1]]), N)
    for (i in 1:N){
      globalAligns1s2 <- pairwiseAlignment(s2, seqaleat[i], substitutionMatrix = M,
                                           gapOpening = openpenalty,
                                           gapExtension = extendedpenalty, scoreOnly = TRUE)

      rep[i] <- globalAligns1s2
    }} else{
      seqaleat <- shuffling(s2, N)
      for (i in 1:N){
        globalAligns1s2 <- pairwiseAlignment(s1, seqaleat[i], substitutionMatrix =M,
                                             gapOpening = openpenalty,
                                             gapExtension = extendedpenalty, scoreOnly = TRUE)
        rep[i] <- globalAligns1s2
      }}

  require("evd")
  ###  Media y desviaciÛn est·ndar de las puntuaciones del vector rep.
  gum_mean = mean(rep)
  gum_sd = sd(rep)

  ##EstimaciÛn  de los par·metros de la distribuciÛn Gumbel.
  gum_lambda =  1.2825 /gum_sd  # scale
  gum_u = gum_mean - 0.45006 * gum_sd  # location

  ##Calculamos la distribucion de Gumble con nuestros parametros
  gumble <- dgumbel(rep, loc = gum_lambda, scale = gum_u)

  ### Estimaci√≥n de K y S'
  K = exp(gum_lambda * gum_u) / (width(s1) * width(s2) )
  S_prima =  gum_lambda * S - log(K *width(seq1) * width(seq2))

  ### Probabilidad
  prob <- sum(rep >= S) / N

  require(Hmisc)

  # 6. RESULTADO GR¡FICO
  require("textplot")

  par(mfrow=c(1,2), oma = c(0, 0, 2, 0))
  hist(rep, main = "Histograma Score", xlab = "Scores",ylab = "Frecuencia de los scores", col = "blue", freq = FALSE)
  curve(dgumbel(x,loc = gum_lambda, scale = gum_u),col="red", from = min(rep), to = max(rep), add = TRUE)
  mtext(paste0("Nombre secuencia 1:",seq1, "\n", "Nombre secuencia 2:", seq2, "\n", "Matriz de substituci√≥n", M, "\n",  "Puntuaci√≥n Gap:",
               openpenalty, extendedpenalty, "\n", "lambda = ",round(gum_lambda,2), "\n", "u=", round(gum_u,2), "\n", "Constante K:", K,
               "\n",  "Score original:", S, "\n", "Score estandarizado (S'):", S_prima, "\n", "P(x >= S)",prob ,"\n"
               "sec_shuffling = ", shuf), side=3, adj = c(1,-1), cex=0.5)
  title(paste("Histograma de los N =", N,"scores"), outer = TRUE, cex.main = 1.5)

## Resultado n˙merico

  txt1 <- "Resumen numÈrico de las puntuaciones obtenidas al hacer shuffling:"
  cat(txt1, fill = TRUE)
  cat(rep("=", nchar(txt1)), sep = "", fill = TRUE)
  print(summary(rep))

  cat("\n\n")
  txt2 <- "Otros resultados numÈricos:"
  cat(txt2, fill = TRUE)
  cat(rep("=", nchar(txt2)), sep = "", fill = TRUE)
  nums <- c(round(gum_lambda,2) , round(gum_u,2), round(K,2), round(S,2), round(S_prima,2), round(prob,2))
  nums <- as.data.frame(nums)
  rownames(nums) <- c("lambda", "u", "K", "S", "S'", "P(x >= S)")
  colnames(nums) <- NULL
  print(nums)
  cat("\n\n")

}

# align -------------------------------------------------------------------

