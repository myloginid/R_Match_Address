rm(list = ls())
options(stringsAsFactors = FALSE) 

check.digits <- function(x){ grepl('^(\\d+)$' , x) }

library(RODBC)
library(RCurl)
library(stringi)
library(stringr)
library(stringdist)
library(sqldf)
library(tm)
library(RTextTools)
library(e1071)
library(plyr)
library(dplyr)

addr.dist = function(x,y, verbose = FALSE) {
# Method name   Description
# s.osa     Optimal string aligment, (restricted Damerau-Levenshtein distance).
# s.lv      Levenshtein distance (as in R's native adist).
# s.dl      Full Damerau-Levenshtein distance.
# s.hamming  Hamming distance (a and b must have same nr of characters).
# s.lcs     Longest common substring distance.
# s.qgram       q-gram distance.
# s.cosine   cosine distance between q-gram profiles
# s.jaccard  Jaccard distance between q-gram profiles
# s.jw      Jaro, or Jaro-Winker distance.
# s.soundex  Distance based on soundex encoding (see below)

    if( exists(x = 'addr.dist', where = -1) == F ) {
        check.digits <- function(x){ grepl('^(\\d+)$' , x) }
    }
    lx =  gsub("[[:punct:]]", "", x )
    ly =  gsub("[[:punct:]]", "", y )
    
    lx =  strsplit(x, split = " ", fixed = TRUE)
    ly =  strsplit(y, split = " ", fixed = TRUE)
    
    if ( verbose == TRUE ) {
        print(lx)
        print(ly)    
    }
    sd.osa  = sd.lv = sd.dl = sd.hamming= sd.lcs = sd.qgram = sd.cosine = sd.jaccard= sd.jw  = sd.soundex= 999.0
    sc.osa  = sc.lv = sc.dl = sc.hamming= sc.lcs = sc.qgram = sc.cosine = sc.jaccard= sc.jw  = sc.soundex= 999.0
    lldx = lldy = llcx = llcy = 0
    
    if ( length(lx[[1]]) > 0  &  length(ly[[1]]) > 0 ) {               
        lxd = unique(  lx[[1]][  check.digits(lx[[1]])   ]  )
        lxc = unique(  lx[[1]][  !check.digits(lx[[1]])  ]  )
        lyd = unique(  ly[[1]][  check.digits(ly[[1]])   ]  )
        lyc = unique(  ly[[1]][  !check.digits(ly[[1]])  ]  )
        sd.osa  = sd.lv = sd.dl = sd.hamming= sd.lcs = sd.qgram = sd.cosine = sd.jaccard = sd.jw  = sd.soundex = 0.0
        llxd = length(lxd)
        llyd = length(lyd)
        llxc = length(lxc)
        llyc = length(lyc)
        
        if ( length(lxd) > 0  &  length(lyd) > 0 ) {        
            if ( verbose == TRUE ) {
                print(lxd)
                print(lyd)                        
            }                    
            for ( k in 1:length(lxd) ) {
                for ( l in 1:length(lyd) ) {
                    if ( verbose == TRUE ) {
                        print(lxd[k])
                        print(lyd[l])                        
                        print(c(sd.osa, sd.dl, sd.dl,sd.hamming,sd.lcs, sd.qgram,sd.cosine,sd.jaccard,sd.jw,sd.soundex))
                    }                    
                    # if ( lxd[k]  == lyd[l] ) { sd.osa          = sd.osa      +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.lv           = sd.lv       +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.dl           = sd.dl       +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.hamming      = sd.hamming  +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.lcs          = sd.lcs      +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.qgram        = sd.qgram    +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.cosine       = sd.cosine   +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.jaccard      = sd.jaccard  +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.jw           = sd.jw       +1 ;  }
                    # if ( lxd[k]  == lyd[l] ) { sd.soundex      = sd.soundex  +1 ;  }
                    if( stringdist( lxd[k], lyd[l] , method = "osa"      ) == 0     )   { sd.osa     = sd.osa     + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "lv"       ) == 0     )   { sd.lv      = sd.lv      + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "dl"       ) == 0     )   { sd.dl      = sd.dl      + 1 ;  } 
                    # if( stringdist( lxd[k], lyd[l] , method = "hamming"  )          )   { sd.hamming = sd.hamming + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "lcs"      ) <  3     )   { sd.lcs     = sd.lcs     + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "qgram"    ) <  3     )   { sd.qgram   = sd.qgram   + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "cosine"   ) <  0.100 )   { sd.cosine  = sd.cosine  + 1 ;  } 
                    # if( stringdist( lxd[k], lyd[l] , method = "jaccard"  )          )   { sd.jaccard = sd.jaccard + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "jw"       ) <  0.100 )   { sd.jw      = sd.jw      + 1 ;  } 
                    if( stringdist( lxd[k], lyd[l] , method = "soundex"  ) == 0     )   { sd.soundex = sd.soundex + 1 ;  }
                }
            }
        }
        sd.osa     = sd.osa      / max( length(lxd) , length(lyd)  )
        sd.lv      = sd.lv       / max( length(lxd) , length(lyd)  )
        sd.dl      = sd.dl       / max( length(lxd) , length(lyd)  )
        sd.hamming = sd.hamming  / max( length(lxd) , length(lyd)  )
        sd.lcs     = sd.lcs      / max( length(lxd) , length(lyd)  )
        sd.qgram   = sd.qgram    / max( length(lxd) , length(lyd)  )
        sd.cosine  = sd.cosine   / max( length(lxd) , length(lyd)  )
        sd.jaccard = sd.jaccard  / max( length(lxd) , length(lyd)  )
        sd.jw      = sd.jw       / max( length(lxd) , length(lyd)  )
        sd.soundex = sd.soundex  / max( length(lxd) , length(lyd)  )

        sc.osa  = sc.lv = sc.dl = sc.hamming= sc.lcs = sc.qgram = sc.cosine = sc.jaccard= sc.jw  = sc.soundex= 0.0
        if ( length(lxc) > 0  &  length(lyc) > 0 ) {        
            if ( verbose == TRUE ) {
                print(lxc)
                print(lyc)                        
            }                    
            for ( k in 1:length(lxc) ) {
                for ( l in 1:length(lyc) ) {
                    if( stringdist( lxc[k], lyc[l] , method = "osa"      ) == 0     )   { sc.osa     = sc.osa     + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "lv"       ) == 0     )   { sc.lv      = sc.lv      + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "dl"       ) == 0     )   { sc.dl      = sc.dl      + 1 ;  } 
                    # if( stringdist( lxc[k], lyc[l] , method = "hamming"  )          )   { sc.hamming = sc.hamming + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "lcs"      ) <  3     )   { sc.lcs     = sc.lcs     + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "qgram"    ) <  3     )   { sc.qgram   = sc.qgram   + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "cosine"   ) <  0.100 )   { sc.cosine  = sc.cosine  + 1 ;  } 
                    # if( stringdist( lxc[k], lyc[l] , method = "jaccard"  )          )   { sc.jaccard = sc.jaccard + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "jw"       ) <  0.100 )   { sc.jw      = sc.jw      + 1 ;  } 
                    if( stringdist( lxc[k], lyc[l] , method = "soundex"  ) == 0     )   { sc.soundex = sc.soundex + 1 ;  }
                    if ( verbose == TRUE ) {
                        print(lxc[k])
                        print(lyc[l])                        
                        print(c(sc.osa, sc.dl, sc.dl,sc.hamming,sc.lcs, sc.qgram,sc.cosine,sc.jaccard,sc.jw,sc.soundex))                    
                    }
                }
            }
        }
        sc.osa     = sc.osa      / max( length(lxc) , length(lyc)  )
        sc.lv      = sc.lv       / max( length(lxc) , length(lyc)  )
        sc.dl      = sc.dl       / max( length(lxc) , length(lyc)  )
        # sc.hamming = sc.hamming  / max( length(lxc) , length(lyc)  )
        sc.lcs     = sc.lcs      / max( length(lxc) , length(lyc)  )
        sc.qgram   = sc.qgram    / max( length(lxc) , length(lyc)  )
        sc.cosine  = sc.cosine   / max( length(lxc) , length(lyc)  )
        # sc.jaccard = sc.jaccard  / max( length(lxc) , length(lyc)  )
        sc.jw      = sc.jw       / max( length(lxc) , length(lyc)  )
        sc.soundex = sc.soundex  / max( length(lxc) , length(lyc)  )
    }
    return ( cbind( llxd , llyd , llxc , llyc
                ,sd.osa, sd.dl, sd.dl,sd.hamming,sd.lcs, sd.qgram,sd.cosine,sd.jaccard,sd.jw,sd.soundex
                ,sc.osa, sc.dl, sc.dl,sc.hamming,sc.lcs, sc.qgram,sc.cosine,sc.jaccard,sc.jw,sc.soundex) )
}

# Some Cleansing, Keep if you need it
# amerged has 2 addresses  - vaddr and caddr. Code below will cleanse, match it and apply the rules to mark them as matched or unmatched


amerged$vaddr = paste(amerged$name_line_2, amerged$name_line_3, amerged$addr_1, amerged$addr_2, amerged$pstl_cd, sep=" ")
amerged$vaddr = tolower(amerged$vaddr)
amerged$vaddr = str_replace(gsub("\\s+", " ", str_trim(amerged$vaddr)), "B", "b")
amerged$vaddr = gsub('(\\s+|^)([0]+)(\\d+)(\\s+|$)','\\1\\3\\4', amerged$vaddr )
amerged$vaddr = gsub( '\\b(\\w+)\\s+\\1\\b','\\1', amerged$vaddr )
amerged$vaddr = gsub( '(singapore)','', amerged$vaddr )
amerged$vaddr = gsub( '(blk)','', amerged$vaddr )

amerged$caddr = paste(amerged$blk, amerged$lvl, amerged$unit, amerged$strtname1, amerged$strtname2, amerged$strtname3, amerged$strtname4, amerged$postcode, sep=" ")
amerged$caddr = tolower(amerged$caddr)
amerged$caddr = str_replace(gsub("\\s+", " ", str_trim(amerged$caddr)), "B", "b")
amerged$caddr = gsub('(\\s+|^)([0]+)(\\d+)(\\s+|$)','\\1\\3\\4', amerged$caddr )
amerged$caddr = gsub( '\\b(\\w+)\\s+\\1\\b','\\1', amerged$caddr )
amerged$caddr = gsub( '999999','', amerged$caddr )

amerged.out = t(mapply(function(x,y) addr.dist(x,y,FALSE), amerged$caddr, amerged$vaddr))

# test = amerged[sample(1:NROW(amerged),round(NROW(amerged)/1000)),]
# mapply(function(x,y) addr.dist(x,y,FALSE), test$caddr, test$vaddr)

row.names(amerged.out) = NULL
colnames(amerged.out) = c(  'ldx', 'ldy', 'lcx', 'lcy'
                            ,'sd.osa', 'sd.dl', 'sd.dl','sd.hamming','sd.lcs', 'sd.qgram','sd.cosine','sd.jaccard','sd.jw','sd.soundex'
                            ,'sc.osa', 'sc.dl', 'sc.dl','sc.hamming','sc.lcs', 'sc.qgram','sc.cosine','sc.jaccard','sc.jw','sc.soundex'
                        )
amerged = cbind(amerged, amerged.out)

amerged$Result = ifelse(  amerged$sd.osa >= 1 & amerged$sc.cosine >= 0.15 , 'Match_R1' , 'UN')
amerged$Result = ifelse(  amerged$Result == '' & amerged$sd.osa >= 1 & amerged$sc.cosine < 0.15 & amerged$ldx>=3 & amerged$ldy>=3  , 'Match_R2' , amerged$Result)
amerged$Result = ifelse(  amerged$Result == '' & amerged$sd.osa == 0 & amerged$sc.cosine >= 0.65 & amerged$ldx==0 & amerged$ldy==0  , 'Match_R3' , amerged$Result)

# store the result where you want
