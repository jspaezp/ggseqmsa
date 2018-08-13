context("plotting")

test_that("multiplication works", {

    fakeseqences_char <- c('PEPTIDES','MYPEPTID','MYPEPTLD',
                           'TIDALESS','PRTMICSR')

    fakesequences_aastrng <- Biostrings::AAStringSet(
        c('PEPTIDES','MYPEPTID','MYPEPTLD',
          'TIDALESS','PRTMICSR'))


    fakeseqences_file <- tempfile()
    Biostrings::writeXStringSet(fakesequences_aastrng, fakeseqences_file)

    cluster.sequences(fakeseqences_char)
    cluster.sequences(fakesequences_aastrng)
    cluster.sequences(fakeseqences_file)

    g<-draw.msa.dendrogam(fakeseqences_char, cutoff = 0.5)
    g<-draw.msa.dendrogam(fakesequences_aastrng, cutoff = 0.5)
    g<-draw.msa.dendrogam(fakeseqences_file, cutoff = 0.5)
    g<-draw.msa.dendrogam(fakeseqences_file, cutoff = 0.5)

})
