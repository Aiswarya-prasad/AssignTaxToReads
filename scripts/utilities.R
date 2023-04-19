make_theme <- function(theme_name=theme_classic() ,max_colors=0, palettefill="Pastel1", palettecolor="Dark2", modify_guide = T,
                        setFill=TRUE, setCol=TRUE,
                        guide_nrow=2, guide_nrow_byrow=TRUE, leg_pos="top", leg_size=12,
                        axis_x_title = 12, axis_y_title = 12,
                        x_angle=0 ,x_vj=0, x_hj=0, x_size=12,
                        y_angle=0 ,y_vj=0, y_hj=0, y_size=12){
  n_11 = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  n_12 = c("Paired", "Set3")
  n_8 = c("Accent", "Dark2", "Pastel2", "Set2")
  if (palettefill %in% n_12) {
    n_f = 12
  } else {
    if (palettefill %in% n_11) {
      n_f = 11
    } else {
      if (palettefill %in% n_8) {
        n_f  = 8
      } else {
        n_f = 9
      }
    }
  }
  if (palettecolor %in% n_12) {
    n_c = 12
  } else {
    if (palettecolor %in% n_11) {
      n_c = 11
    } else {
      if (palettecolor %in% n_8) {
        n_c  = 8
      } else {
        n_c = 9
      }
    }
  }
  getFill = colorRampPalette(brewer.pal(n_f, palettefill))
  getColor = colorRampPalette(brewer.pal(n_c, palettecolor))
  theme_params <- theme(axis.text.x = element_text(angle = x_angle,
    vjust = x_vj, hjust=x_hj,
    size = x_size),
    axis.text.y = element_text(angle = y_angle,
      vjust = y_vj, hjust=y_hj,
      size = y_size),
      axis.title.x = element_text(size=axis_x_title),
      axis.title.y = element_text(size=axis_y_title),
      legend.position=leg_pos,
      legend.text = element_text(size=leg_size)
    )
  if (modify_guide == T) {
    guide_params <- guides(fill = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  ),
                          col = guide_legend(
                                    nrow=guide_nrow,
                                    byrow=guide_nrow_byrow
                                  )
                    )
  my_theme <- list(
                theme_name,
                theme_params,
                guide_params
              )
  } else {
    my_theme <- list(
                  theme_name,
                  theme_params
                )
  }
  if(setFill) {
    if (n_f < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_fill_manual(values = getFill(max_colors), na.value="grey")
                  )
    } else {
      my_theme <- list(
                    my_theme,
                    scale_fill_brewer(palette=palettefill, na.value="grey")
                  )
    }
  }
  if(setCol) {
    if (n_c < max_colors) {
      my_theme <- list(
                    my_theme,
                    scale_color_manual(values = getColor(max_colors), na.value="grey")
                  )
    } else {
      my_theme <- list(
                    my_theme,
                    scale_color_brewer(palette=palettecolor, na.value="grey")
                  )
    }
  }
  return(my_theme)
}

remove_extension <- function(x, extension) {
  strsplit(x, extension)[[1]][[1]]
}

get_reads_path <- function(ID) {
  return(paste0(raw_reads_renamed_path, ID, "_reads.fastq.gz"))
}
get_number_of_reads <- function(file) {
  commmand <- paste0("zcat ",file,"| grep -e'^+$' | wc -l")
  number <- system(commmand, intern = TRUE)
  return(as.numeric(number))

}
get_number_assigned <- function(rank_name) {
  if (rank_name == "Total") {
    return(dim(taxonomy_df)[[1]])
  } else {
    num_ASVs <- dim(taxonomy_df[rank_name])[[1]] - sum(is.na(taxonomy_df[rank_name]))
    return(num_ASVs)
  }
}
genusColors <- c("Bombilactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[1],
                    "Lactobacillus" = head(colorRampPalette(c(brewer.pal(11, "Spectral")[1], "#FFFFFF"))(10), -1)[4],
                    "Bifidobacterium" = brewer.pal(11, "Spectral")[3],
                    "Gilliamella" = brewer.pal(11, "Spectral")[11],
                    "Frischella" = brewer.pal(11, "Spectral")[8],
                    "Bartonella" = brewer.pal(11, "Spectral")[7],
                    "Snodgrassella" = brewer.pal(11, "Spectral")[10],
                    "Apibacter" = brewer.pal(11, "Spectral")[4],
                    "Commensalibacter" = brewer.pal(11, "Spectral")[6],
                    "Bombella" = brewer.pal(11, "Spectral")[5],
                    "Apilactobacillus" = brewer.pal(11, "Spectral")[9],
                    "Dysgonomonas" = brewer.pal(11, "Spectral")[2],
                    "Spiroplasma" = brewer.pal(8, "Set1")[8],
                    "WRHT01" = brewer.pal(8, "Dark2")[3],
                    "Pectinatus" = brewer.pal(8, "Dark2")[1],
                    "Enterobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[1],
                    "Zymobacter" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[2],
                    "Entomomonas"= head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[4],
                    "Saezia" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[6],
                    "Parolsenella" = head(colorRampPalette(c(brewer.pal(11, "BrBG")[2], "#FFFFFF"))(10), -1)[8]
)
extend_colors <- function(names_vec, colors_vec, greys = T, pal = "Pastel1"){
  final_list <- c()
  if (greys) {
     for (a_name in names_vec) {
      if (a_name %in% names(colors_vec)) {
        final_list[a_name] = colors_vec[a_name]
      } else {
        final_list[a_name] = "grey"
      }
    }
  } else {
    i = 1
    num_new_cols = length(names_vec[which(!(names_vec %in% names(colors_vec)))])
    for (a_name in names_vec) {
      if (a_name %in% names(colors_vec)) {
        final_list[a_name] = colors_vec[a_name]
      } else {
        if (num_new_cols > 9) {
         final_list[a_name] = colorRampPalette(brewer.pal(9, pal))(num_new_cols)[i] 
        } else {
          final_list[a_name] = brewer.pal(num_new_cols, pal)[i] 
        }
        i = i + 1
      }
    }
  }
  return(final_list) 
}
