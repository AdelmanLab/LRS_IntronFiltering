#Filtering pipeline used to generate the working intron annotations for uninduced and induced cell conditions related to Reimer et al., 2020   
  #Input File: 
    #Intron annotation derived from active annotated transcription start sites with combined replicate LRS counts per intron. 
    #The uninduced (Final_untreated_spliced, Final_untreated_unspliced) and induced (Final_DMSO_spliced, Final_DMSO_unspliced) LRS counts per intron are in different columns
    #The annotated intron number for minus strand annotations was corrected with respect to the gene start i.e intron 0 corresponds to the intron closest to the TSS and not the terminal intron. 

intron_start <- read.delim("intron_splicing_counts_for_PROseq_classification_uniqueIDs_Condition_Counts_Reported_intron_num_corrected.txt", as.is = T)

#Remove annotations with an intron length greather than 10kb 
intron_start$Intron_Size <- abs(intron_start$i_start - intron_start$i_end)
intron_lengthfilter <- subset(intron_start, intron_start$Intron_Size <= 10000)

#Remove annotations with 0 spliced read counts per condition 
intron_norm_MEL <- subset(intron_lengthfilter, intron_lengthfilter$Final_Untreated_spliced > 0)
intron_norm_DMSO <- subset(intron_lengthfilter, intron_lengthfilter$Final_DMSO_spliced > 0)

write.table(intron_norm_MEL, "introns_MEL_spliced>0.txt", sep="\t", row.names = FALSE, quote = FALSE)
write.table(intron_norm_DMSO, "introns_DMSO_spliced>0.txt", sep="\t", row.names = FALSE, quote = FALSE)

#Rank order introns by strand, Chr, 5'SS Coordinate, 3'SS Coordinate, and Intron Number (corrected column) 
intron_MEL_spliced_sorted <- intron_norm_MEL[order(intron_norm_MEL$i_strand, intron_norm_MEL$i_chr, intron_norm_MEL$FiveSS, intron_norm_MEL$ThreeSS, intron_norm_MEL$intron_num_corrected) , ]
intron_DMSO_spliced_sorted <- intron_norm_DMSO[order(intron_norm_DMSO$i_strand, intron_norm_DMSO$i_chr, intron_norm_DMSO$FiveSS, intron_norm_DMSO$ThreeSS, intron_norm_DMSO$intron_num_corrected) , ]

#Filtering loop to do the following...
  # 1. For annotation duplicates, the annotation with the smallest intron number was retained
  # 2. If two annotations share a 5’ or 3’SS, the intron with the most spliced counts was retained. 

#To generate the uninduced intron list...
#loop variables 
n=0                         # to count how many times the loop was repeated 
nrowremoved = 1             # Condition term for the while loop. Keep subsetting until there are no more introns (rows) to remove... 

loop_temp <- intron_MEL_spliced_sorted
loop_temp_flagged <- read.delim("intron_annotations_flagged.txt", as.is = T) #Empty starting file to retain flagged annotations 

#Filtering Loop: 
while ( nrowremoved !=0 ) {
  
  # For the first loop iteration, make two new column or 
  if ( n == 0 ){
    loop_temp$flag <- 0               #Flag Terms. 0 = keep, 1+ = remove
    loop_temp_flagged$flag <- 0           
  
  #Reset the loop starting file by removing flagged introns from the previous iteration  
  } else {
    loop_temp <- loop_temp_final      #loop_temp_final is made at the end of the loop and is a subset of all annotations with a flag of 0
    rownames(loop_temp) <- NULL
  }
  
  #For each row in the document...
  for (i in 1:dim(loop_temp)[1]) {
    
    #Condition: If intron[i] and intron[i+1] have the same 5'SS and 3'SS, Add a +1 Flag to the annotation with a larger intron number
    
    if( isTRUE( (loop_temp$FiveSS[i] == loop_temp$FiveSS[i+1]) & (loop_temp$ThreeSS[i] == loop_temp$ThreeSS[i+1]) )) {        
      
      if(isTRUE( loop_temp$intron_num[i] > loop_temp$intron_num[i+1]) ) {     
        loop_temp$flag[i] = loop_temp$flag[i] + 1
      } else {            
        loop_temp$flag[i+1] = loop_temp$flag[i+1] + 1
      }
    
    #Condition: If intron[i] and intron[i+1] have the same 5'SS or 3'SS, Add a +1 Flag to the annotation with fewer spliced LRS read counts 
      
    } else if ( isTRUE(((loop_temp$FiveSS[i] == loop_temp$FiveSS[i+1]) &  (loop_temp$ThreeSS[i] != loop_temp$ThreeSS[i+1])) | ((loop_temp$FiveSS[i] != loop_temp$FiveSS[i+1]) &  (loop_temp$ThreeSS[i] == loop_temp$ThreeSS[i+1]))) ) {            
      
      if(isTRUE( loop_temp$Final_Untreated_spliced[i] < loop_temp$Final_Untreated_spliced[i+1]) ) {     
        loop_temp$flag[i] = loop_temp$flag[i] + 1
      } else {            
        loop_temp$flag[i+1] = loop_temp$flag[i+1] + 1
      }
    } else {
    } 
  }
  
  #Condition for the while loop. 
  loop_temp_final <- subset(loop_temp, loop_temp$flag == 0 )                                #Subset transcripts with a flag=0 
  loop_temp_flagged <- rbind(loop_temp_flagged, subset(loop_temp, loop_temp$flag > 0 ))     #Add flagged intron annotations to the flagged output file 
  nrowremoved <- nrow(loop_temp) - nrow(loop_temp_final)                                    #While loop variable. If 0, the while loop will end meaning that all flagged annotations have been removed  
  n=n+1                                                                                     #Track loop iteration
  
}

#Output Files 
write.table(loop_temp_final, "Working_Intron_Annotation_MEL.txt", sep="\t", row.names = FALSE )
write.table(loop_temp_flagged, "Flagged_Intron_Annotation_MEL.txt", sep="\t", row.names = FALSE )

#Repeat code (Lines 28-86) to generate the induced intron list...
n=0                         
nrowremoved = 1             
loop_temp <- intron_DMSO_spliced_sorted
loop_temp_flagged <- read.delim("intron_annotations_flagged.txt", as.is = T) 

#Filtering Loop: 
while ( nrowremoved !=0 ) {
  
  if ( n == 0 ){
    loop_temp$flag <- 0               
    loop_temp_flagged$flag <- 0           

  } else {
    loop_temp <- loop_temp_final    
    rownames(loop_temp) <- NULL
  }
  
  for (i in 1:dim(loop_temp)[1]) {
    if( isTRUE( (loop_temp$FiveSS[i] == loop_temp$FiveSS[i+1]) & (loop_temp$ThreeSS[i] == loop_temp$ThreeSS[i+1]) )) {        
      
      if(isTRUE( loop_temp$intron_num[i] > loop_temp$intron_num[i+1]) ) {     
        loop_temp$flag[i] = loop_temp$flag[i] + 1
      } else {            
        loop_temp$flag[i+1] = loop_temp$flag[i+1] + 1
      }
      
    } else if ( isTRUE(((loop_temp$FiveSS[i] == loop_temp$FiveSS[i+1]) &  (loop_temp$ThreeSS[i] != loop_temp$ThreeSS[i+1])) | ((loop_temp$FiveSS[i] != loop_temp$FiveSS[i+1]) &  (loop_temp$ThreeSS[i] == loop_temp$ThreeSS[i+1]))) ) {            
      
      if(isTRUE( loop_temp$Final_DMSO_spliced[i] < loop_temp$Final_DMSO_spliced[i+1]) ) {     
        loop_temp$flag[i] = loop_temp$flag[i] + 1
      } else {            
        loop_temp$flag[i+1] = loop_temp$flag[i+1] + 1
      }
      
    } else {
    } 
  }

  loop_temp_final <- subset(loop_temp, loop_temp$flag == 0 )                            
  loop_temp_flagged <- rbind(loop_temp_flagged, subset(loop_temp, loop_temp$flag > 0 )) 
  nrowremoved <- nrow(loop_temp) - nrow(loop_temp_final)                                  
  n=n+1                                                                                 

}

write.table(loop_temp_final, "Working_Intron_Annotation_DMSO.txt", sep="\t", row.names = FALSE )
write.table(loop_temp_flagged, "Flagged_Intron_Annotation_DMSO.txt", sep="\t", row.names = FALSE )




