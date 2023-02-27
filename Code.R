library("seqinr")
library("seqinr")
# Reading sequence data into R

dengueseq <- read.fasta("sequence.fasta")
dengueseq
dengue <- read.fasta(file = "sequence.fasta")
dengueseq <- dengue[[1]]

#Length of a DNA sequence
length(dengueseq)


#Base composition of a DNA sequence
table(dengueseq)

#GC Content of DNA
(2240+2770)*100/(3426+2240+2770+2299)
GC(dengueseq)

#DNA words
count(dengueseq, 1)

#As expected, this gives us the number of occurrences of the individual nucleotides. To find the number of occurrences of DNA words that are 2 nucleotides long, we type:
count(dengueseq, 2)

#Note that by default the count() function includes all overlapping DNA words in a sequence. Therefore, for example, the sequence “ATG” is considered to contain two words that are two nucleotides long: “AT” and “TG”.

#If you type help(‘count’), you will see that the result (output) of the function count() is a table object. This means that you can use double square brackets to extract the values of elements from the table. For example, to extract the value of the third element (the number of Gs in the DEN-1 Dengue virus sequence), you can type:
denguetable <- count(dengueseq,1)
denguetable[[3]]

#The command above extracts the third element of the table produced by count(dengueseq,1), which we have stored in the table variable denguetable.

#Alternatively, you can find the value of the element of the table in column “g” by typing:
denguetable[["g"]]



#The following command prints out the first 50 nucleotides of the DEN-1 Dengue virus genome sequence:

dengueseq[1:50]

