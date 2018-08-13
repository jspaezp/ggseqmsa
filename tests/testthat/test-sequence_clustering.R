context("sequence_clustering")

test_that("multiplication works", {
    fakeseqences_char <- c('PEPTIDES','MYPEPTID','MYPEPTLD',
                           'TIDALESS','PRTMICSR')

    fakesequences_aastrng <- Biostrings::AAStringSet(
        c('PEPTIDES','MYPEPTID','MYPEPTLD',
          'TIDALESS','PRTMICSR'))


    fakeseqences_file <- tempfile()
    Biostrings::writeXStringSet(fakesequences_aastrng, fakeseqences_file)
  expect_equal(2 * 2, 4)
})
