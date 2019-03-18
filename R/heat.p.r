
heat.p <- function (df, p.breaks = c(0.001, 0.01, 0.05), 
                    colors = c("darkblue", "blue", "lightblue", "grey"),
                    title = NULL, subtitle = NULL,
                    filter = c(FALSE, FALSE), x.size, y.size, symbol.size,
                    zigzag = c(FALSE, FALSE), abbrv = c(FALSE, FALSE), 
                    margin = c(1, 1), legend = TRUE) 
{
  if (!requireNamespace("ggplot2")) install.packages("ggplot2")
  if (!requireNamespace("grid")) install.packages("grid")
  require(ggplot2)
  require(grid)
  
  df = as.data.frame(df)
  if (!"responses" %in% names(df) | !"variables" %in% names(df) | 
      !"pvalue" %in% names(df) | !"Estimate" %in% names(df) ) 
    stop("df should include 'responses', 'variables', and 'pvalue'")
  if (length(colors) != length(p.breaks) + 1) 
    stop("number of colors does not equal p value intervals")
  
  filter.res = filter[1] 
  filter.var = filter[2]
  y.zigzag = zigzag[1]
  x.zigzag = zigzag[2]
  y.abbrv = abbrv[1]
  x.abbrv = abbrv[2]
  left.margin = margin[1]
  bottom.margin = margin[2]
  
  df = df[df[, "variables"] != "(Intercept)", ]
  br = c(0, p.breaks, 1)
  br = unique(sort(br))
  df$p_value = cut(df$pvalue, br)
  if (filter.res) {
    keep.res = as.character(unique(df[df[, "pvalue"] < br[length(br) - 
                                                            1], "responses"]))
    df = df[df[["responses"]] %in% keep.res, ]
  }
  if (filter.var) {
    keep.var = as.character(unique(df[df[, "pvalue"] < br[length(br) - 
                                                            1], "variables"]))
    df = df[df[["variables"]] %in% keep.var, ]
  }
  if (missing(x.size)) 
    x.size <- NULL
  if (missing(y.size)) 
    y.size <- NULL
  p = ggplot(df, aes(variables, responses))
  p = p + xlab("") + ylab("")
  p = p + geom_tile(aes(fill = p_value))
  p = p + scale_fill_manual(values = colors)
  p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                           hjust = 1, size = x.size))
  p = p + theme(axis.text.y = element_text(size = y.size))
  if (legend) p = p + theme(legend.justification = "top")
  else p = p + theme(legend.position = "none")
  if(x.abbrv) p = p + scale_x_discrete(label=abbreviate)
  if(y.abbrv) p = p + scale_y_discrete(label=abbreviate)
  df.sub = df[df[, "Estimate"] > 0 & df[, "pvalue"] < br[length(br) - 1], ]
  if (missing(symbol.size)) 
    p = p + geom_text(data = df.sub, aes(variables, responses, 
                                         label = "+"), col = "white")
  else p = p + geom_text(data = df.sub, aes(variables, responses, 
                                            label = "+"), col = "white", size = symbol.size)
  p = p + labs(title=title, subtitle=subtitle)
  
  g = ggplotGrob(p)
  g$widths[3] = unit(left.margin, "cm")
  g$heights[7] = unit(bottom.margin, "cm")
  if(y.zigzag){
    g = editGrob(grid.force(g),
               gPath("axis-l", "axis", "axis", "GRID.text"),
               x = unit(c(1,0), "npc"),
               grep = TRUE)

  }
  if(x.zigzag){
    g = editGrob(grid.force(g),
               gPath("axis-b", "axis", "axis", "GRID.text"),
               y = unit(c(1,0), "npc"),
               grep = TRUE)
  }

  grid.newpage()
  grid.draw(g)
  return(g)
}

