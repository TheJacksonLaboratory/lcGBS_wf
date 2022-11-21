library(dplyr)
library(tidyr)
library(parallel)
library(vroom)
library(qtl2convert)
library(furrr)

# load in 4WC genotypes pulled from VCF
load("data/04_pull_4WC_founder_genos_out.RData")
# as of 20221118 genotypes are *exact* position matches

# read in consensus genotypes for PWD and CAST (PWK obviously not a perfect proxy for PWK)
consensus_genotypes <- list.files("data/GM/",pattern = "GM_foundergeno")

# read in allele codes and founder genotype files as templates
GM_allele_codes <- vroom::vroom("data/GM/GM_allelecodes.csv", skip = 3)
example_founder_geno <- vroom::vroom(file = paste0("data/GM/", consensus_genotypes[[1]]), skip = 3)

# clean up genotypes
if(nrow(data.frame(rownames(genos))) == nrow(genos)){
  new_genos <- cbind(data.frame(rownames(genos)),genos)
  rownames(new_genos) <- NULL
}
tidy_4WC_genotypes <- new_genos %>%
  dplyr::filter(count.fields(textConnection(new_genos$rownames.genos.), sep = ";") == 1) %>%
  tidyr::separate(rownames.genos., c("chr","pos","REF","ALT"))

# clean up supplied marker information
marker_info <- supplied_gm_metadata %>%
  dplyr::rename(chr = seqnames,
                pos = start) %>%
  dplyr::select(marker, chr, pos) %>%
  dplyr::mutate(pos = as.character(pos))

# combine 4WC genotypes from VCF and marker information where common
marker_genotypes_4WC <- dplyr::left_join(tidy_4WC_genotypes, marker_info)
head(marker_genotypes_4WC)

# Generate allele codes
allele_codes_4WC <- marker_genotypes_4WC %>%
  dplyr::select(marker, chr, REF, ALT) %>%
  dplyr::rename(A = REF,
                B = ALT)

# Take standing 4WC allele codes and recode to match what is present in the GM allele codes file
final_allele_codes_4WC <- allele_codes_4WC %>%
  dplyr::rename(A_4WC = A,
                B_4WC = B) %>%
  # Join to GM allele codes
  dplyr::left_join(., GM_allele_codes) %>%
  # Adjust which nucleotide = (A|B) in allele codes
  dplyr::mutate(A_4WC = dplyr::case_when(is.na(A) ~ A_4WC,
                                         A_4WC == A ~ A_4WC,
                                         A_4WC != A ~ A),
                B_4WC = dplyr::case_when(is.na(B) ~ B_4WC,
                                         B_4WC == B ~ B_4WC,
                                         B_4WC != B ~ B)) %>%
  dplyr::select(-A, -B) %>%
  dplyr::rename(A = A_4WC, 
                B = B_4WC) %>%
  dplyr::distinct()


if(!file.exists("data/4WC_founder_genos/")){
  dir.create(path = "data/4WC_founder_genos")
}

# Generate 4WC founder genotypes for each chromosome
# > colnames(example_founder_geno)
# [1] "marker" "A"      "B"      "C"      "D"      "E"      "F"      "G"      "H"

# Nest VCF founder genotypes by chromosome
chr_nested_marker_genotypes_4WC <- marker_genotypes_4WC %>%
  dplyr::select(-pos) %>%
  dplyr::distinct() %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()

# Nest marker allele codes by chromosome
chr_nested_allele_codes_4WC <- final_allele_codes_4WC %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()

genos <- chr_nested_marker_genotypes_4WC$data[[2]]
alleleCodes <- chr_nested_allele_codes_4WC$data[[2]]
chromosome <- chr_nested_marker_genotypes_4WC$chr[[2]]
## Function to convert VCF genotype codes to founder genotype files
assign4WCFounderGenos <- function(genos, alleleCodes, chromosome){
  # Recode VCF genotypes to letter codes
  nucleotideRecoded <- genos %>%
    dplyr::select(marker, everything()) %>%
    dplyr::mutate(CAST = dplyr::case_when(CAST == "1/1" ~ ALT,
                                          CAST == "0/0" ~ REF,
                                          CAST == "0/1" ~ "HET",
                                          CAST == "./." ~ "-"),
                  GOR = dplyr::case_when(GOR == "1/1" ~ ALT,
                                         GOR == "0/0" ~ REF,
                                         GOR == "0/1" ~ "HET",
                                         GOR == "./." ~ "-"),
                  POHN = dplyr::case_when(POHN == "1/1" ~ ALT,
                                          POHN == "0/0" ~ REF,
                                          POHN == "0/1" ~ "HET",
                                          POHN == "./." ~ "-"),
                  PWD = dplyr::case_when(PWD == "1/1" ~ ALT,
                                         PWD == "0/0" ~ REF,
                                         PWD == "0/1" ~ "HET",
                                         PWD == "./." ~ "-")) %>%
    # Remove sites that segregate within founder strains (can discuss this more....)
    dplyr::mutate(het = if_else(CAST == "HET" | GOR == "HET" | POHN == "HET" | PWD == "HET", true = "HET", false = "NO_HET")) %>%
    dplyr::filter(het != "HET") %>%
    dplyr::select(-het, -REF, -ALT)
  
  # Use allele codes for all genotypes
  alleleRecoded <- nucleotideRecoded %>%
    dplyr::inner_join(.,alleleCodes) %>%
    dplyr::distinct() %>%
    dplyr::mutate(CAST = dplyr::if_else(CAST == A, true = "A", false = "B"),
                  GOR = dplyr::if_else(GOR == A, true = "A", false = "B"),
                  POHN = dplyr::if_else(POHN == A, true = "A", false = "B"),
                  PWD = dplyr::if_else(PWD == A, true = "A", false = "B")) %>%
    dplyr::select(-A, -B)
    
  # Substitute letter codes as column names like GM example file
  colnames(alleleRecoded)[colnames(alleleRecoded) == "CAST"] <- "F"
  colnames(alleleRecoded)[colnames(alleleRecoded) == "POHN"] <- "J"
  colnames(alleleRecoded)[colnames(alleleRecoded) == "GOR"] <- "K"
  colnames(alleleRecoded)[colnames(alleleRecoded) == "PWD"] <- "L"
  
  # Re-order strains in alphabetical order and remove duplicate or NA markers
  alleleRecodedStrainLetters <- alleleRecoded[!duplicated(alleleRecoded$marker),] %>%
    dplyr::select(marker, F, J, K, L) %>%
    dplyr::filter(!is.na(marker))

  # Write csv in qtl2 format
  qtl2convert::write2csv(df = alleleRecodedStrainLetters, 
                         filename = paste0("data/4WC_founder_genos/4WC_foundergeno",chromosome,".csv"),
                         comment = paste0("Founder genotypes (A/B/-) for 4WC; chr ", chromosome),
                         overwrite = T)
}
purrr::pmap(.l = list(chr_nested_marker_genotypes_4WC$data,
                      chr_nested_allele_codes_4WC$data,
                      chr_nested_marker_genotypes_4WC$chr), 
            .f = assign4WCFounderGenos)

