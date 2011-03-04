library(yasomi)
library(RUnit)
testsuite.somgrid <- defineTestSuite("somgrid",dirs="unit",
                                     testFileRegexp="^runitsomgrid.*\\.[rR]$")
testresult.somgrid <- runTestSuite(testsuite.somgrid)
printTextProtocol(testresult.somgrid)
