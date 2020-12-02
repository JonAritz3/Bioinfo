#' Aleatorizacion de una secuencia.
#'
#' La funcion \code{generatessequences} genera un numero determinado (por ejemplo, 1000) de 
#' secuencias aleatorias de nucleotidos (si la secuencia es ADN) o de aminoacidos (si la secuencia es una
#' proteina) de una longitud determinada mediante un modelo multinomial.
#' 
#' @param inputsequence objeto tipo cadena de carácteres con la secuencia a aleatorizar.
#' @param X obtejo numerico. Numero de veces que se desea aleatorizar la secuencia.
#'
#' @return  \code{generateSeqsWithMultinomialModel} devuelve un vector de longitud X con elementos de tipo
#' caracter, que son las secuencias ordenadas al azar.
#'
#' @examples
#' 
#' generateSeqsWithMultinomialModel("PAWHEAE",100)
#' generateSeqsWithMultinomialModel("ATTGGCACT",10)
#' 
#' @references \url{https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter4.html}
#'
#' @export
#' 

shuffling <- function(inputsequence, X){
  require("seqinr")
  inputsequencevector <- s2c(inputsequence)
  mylength <- length(inputsequencevector)
  mytable <- table(inputsequencevector)
  letters <- rownames(mytable)
  numletters <- length(letters)
  probabilities <- numeric()
  for (i in 1:numletters)
  {
    letter <- letters[i]
    count <- mytable[[i]]
    probabilities[i] <- count/mylength
  }
  seqs <- numeric(X)
  for (j in 1:X)
  {
    seq <- sample(letters, mylength, rep=TRUE, prob=probabilities)
    seq <- c2s(seq)
    seqs[j] <- seq
  }
  return(seqs)
}

