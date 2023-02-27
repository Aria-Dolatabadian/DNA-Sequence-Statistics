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

#A sliding window analysis of GC content
GC(dengueseq[1:2000])      # Calculate the GC content of nucleotides 1-2000 of the Dengue genome

GC(dengueseq[2001:4000])   # Calculate the GC content of nucleotides 2001-4000 of the Dengue genome

GC(dengueseq[4001:6000])   # Calculate the GC content of nucleotides 4001-6000 of the Dengue genome

GC(dengueseq[6001:8000])   # Calculate the GC content of nucleotides 6001-8000 of the Dengue genome

GC(dengueseq[8001:10000])  # Calculate the GC content of nucleotides 8001-10000 of the Dengue genome

GC(dengueseq[10001:10735]) # Calculate the GC content of nucleotides 10001-10735 of the Dengue genome


#Or

starts <- seq(1, length(dengueseq)-2000, by = 2000)
starts

n <- length(starts)    # Find the length of the vector "starts"
for (i in 1:n) {
        chunk <- dengueseq[starts[i]:(starts[i]+1999)]
        chunkGC <- GC(chunk)
        print (chunkGC)
     }


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

#It is common to use the data generated from a sliding window analysis to create a sliding window plot of GC content. 
#To create a sliding window plot of GC content, you plot the local GC content in each window of the genome, versus the nucleotide position of the start of each window. 
#We can create a sliding window plot of GC content by typing:

starts <- seq(1, length(dengueseq)-2000, by = 2000)
n <- length(starts)    # Find the length of the vector "starts"
chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
for (i in 1:n) {
        chunk <- dengueseq[starts[i]:(starts[i]+1999)]
        chunkGC <- GC(chunk)
        print(chunkGC)
        chunkGCs[i] <- chunkGC
     }
plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")


x11()


slidingwindowplot <- function(windowsize, inputseq)
{
   starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
   n <- length(starts)    # Find the length of the vector "starts"
   chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
   for (i in 1:n) {
        chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
        chunkGC <- GC(chunk)
        print(chunkGC)
        chunkGCs[i] <- chunkGC
   }
   plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}

slidingwindowplot(3000, dengueseq)

x11()

slidingwindowplot(300, dengueseq)





count(dengueseq, 1) # Get the number of occurrences of 1-nucleotide DNA words

2770/(3426+2240+2770+2299) # Get fG

2240/(3426+2240+2770+2299) # Get fC

count(dengueseq, 2) # Get the number of occurrences of 2-nucleotide DNA words
  
500/(1108+720+890+708+901+523+261+555+976+500+787+507+440+497+832+529) # Get fGC

0.04658096/(0.2580345*0.2086633) # Get rho(GC)
