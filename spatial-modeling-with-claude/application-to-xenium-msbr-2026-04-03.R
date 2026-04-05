# ---- set wd ----
setwd("~/Desktop/spatial-modeling-with-claude")

# ---- set up ----
require(Seurat)
require(SeuratObject)
library(sf)
library(ggplot2)
library(patchwork)
source("fit_spatial_mask.R")
source("estimate_concentration_field.R")
set.seed(2024)

# ---- shared theme ----
theme_demo <- function(bs=9.5)
  theme_minimal(base_size=bs) +
  theme(plot.title=element_text(face="bold",size=bs+0.5),
        plot.subtitle=element_text(size=bs-1,color="grey40",margin=margin(b=4)),
        panel.grid=element_line(color="grey94"),
        legend.key.width=unit(0.35,"cm"),
        legend.title=element_text(size=bs-1),
        legend.text=element_text(size=bs-1.5),
        axis.text=element_text(size=bs-2))

plot_field <- function(result, transcripts=NULL, title="", subtitle="",
                       palette="magma", log_scale=FALSE,
                       show_pts=TRUE, show_contours=TRUE, n_contours=6L,
                       pt_size=0.35, pt_alpha=0.45, pt_color="#00e5ff",
                       fill_label="Conc.", symmetric=FALSE) {
  df       <- field_to_df(result, inside_only=TRUE)
  fill_col <- if (log_scale) { df$fv <- log1p(df$field); "fv" } else "field"
  fill_lbl <- if (log_scale) paste0("log\u2081\u208a(",fill_label,")") else fill_label
  p <- ggplot(df, aes(x=x,y=y,fill=.data[[fill_col]])) +
    geom_raster(interpolate=TRUE) + coord_equal() +
    labs(title=title, subtitle=subtitle, x="X", y="Y") + theme_demo()
  if (symmetric) {
    lim <- max(abs(df[[fill_col]]),na.rm=TRUE)
    p   <- p + scale_fill_distiller(palette="RdBu",limits=c(-lim,lim),
                                    name=fill_lbl,na.value="transparent")
  } else {
    p <- p + scale_fill_viridis_c(option=palette,name=fill_lbl,na.value="transparent")
  }
  if (show_contours && any(!is.na(df$field)))
    p <- p + geom_contour(aes(z=field),color="white",
                          alpha=0.35,linewidth=0.25,bins=n_contours)
  if (show_pts && !is.null(transcripts))
    p <- p + geom_point(data=transcripts,aes(x=x,y=y),inherit.aes=FALSE,
                        color=pt_color,size=pt_size,alpha=pt_alpha)
  p
}

sample_in_mask <- function(mask_poly, n, max_tries=20) {
  bb <- sf::st_bbox(mask_poly); mu <- sf::st_union(mask_poly)
  pts <- data.frame(x=numeric(0),y=numeric(0))
  for (i in seq_len(max_tries)) {
    if (nrow(pts)>=n) break
    cnd <- data.frame(x=runif(n*8,bb["xmin"],bb["xmax"]),
                      y=runif(n*8,bb["ymin"],bb["ymax"]))
    ok  <- as.logical(sf::st_within(
      sf::st_as_sf(cnd,coords=c("x","y")),mu,sparse=FALSE)[,1L])
    pts <- rbind(pts,cnd[ok,])
  }
  pts[seq_len(min(n,nrow(pts))),]
}

# ---- Shared helpers ----------------------------------------------------------

plot_mask <- function(mask, coords, title = "", subtitle = "",
                      fill = "#4c9be8", border = "#1a5fa8",
                      pt_col = "#c0392b", pt_size = 0.9, pt_alpha = 0.5) {
  # geom_sf is required for correct hole rendering.
  # geom_polygon draws a line between consecutive rings in the same group,
  # producing diagonal artefacts across hole interiors.
  mask_sf <- sf::st_sf(geometry = mask)
  ggplot() +
    geom_sf(data = mask_sf,
            fill = fill, alpha = 0.2, color = border, linewidth = 0.55) +
    geom_point(data = coords, aes(x = x, y = y),
               color = pt_col, size = pt_size, alpha = pt_alpha) +
    coord_sf(datum = NA) +   # datum=NA: suppress geographic graticule
    labs(title = title, subtitle = subtitle, x = "X", y = "Y") +
    theme_minimal(base_size = 9.5) +
    theme(plot.title    = element_text(face = "bold", size = 9.5),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          panel.grid    = element_line(color = "grey93"))
}

# ---- apply to real world data ----
# load xenium object
load("~/Large Files/Pneumonectomy_Project_transfer/MSBR_polishing_2025-11-01/outputs/PPLR.polished.with.territory.and.neighborhood.2025-11-17.Robj")

# subset out just one TMA puck
PPLR.polished <- UpdateSeuratObject(PPLR.polished)
Idents(PPLR.polished) <- PPLR.polished$TMA
table(Idents(PPLR.polished))
temp <- subset(PPLR.polished,
               idents = 'PNX7_M1_Sp')
temp

# get centroid xy coordinates
tissue_pts <- data.frame(x = temp$x,
                      y = temp$y)
# make tissue mask
tissue_mask <- fit_spatial_mask(tissue_pts, 
                                method="raster",
                                raster_resolution=256L, 
                                raster_sigma=1, 
                                raster_threshold=0.75, 
                                verbose=TRUE,
                                buffer_dist = 3,
                                smooth_mask = T,
                                smooth_tolerance = 3)
# plot mask
plot_mask(tissue_mask, 
          tissue_pts,
          title    = "PNX7_M1_Sp",
          subtitle = "raster_resolution = 300, raster_sigma = 1, raster_threshold = 0.75, smooth = T, smooth_tolerance = 3, buffer = 3")

# get transcript xy coordinates
molecules <- temp[['PNX7_M1_Sp']][['molecules']]
tc_Areg <- molecules$Areg@coords
tc_Vegfa <- molecules$Vegfa@coords

# inspect
tissue_sf <- sf::st_sf(geometry = tissue_mask)
print(ggplot() +
        geom_sf(data=tissue_sf, fill="grey92", color="grey60", linewidth=0.4) +
        geom_point(data=tc_Areg,aes(x=x,y=y,color="Areg"),alpha=0.6,size=0.1) +
        geom_point(data=tc_Vegfa,aes(x=x,y=y,color="Vegfa"),alpha=0.5,size=0.1) +
        scale_color_manual(values=c("Areg"="#e74c3c","Vegfa"="#2980b9"),name=NULL) +
        #scale_size_continuous(range=c(0.3,2.5),name="log(UMI+1)") +
        coord_sf(datum=NA) + theme_demo() +
        labs(title="PNX7_M1_Sp",
             subtitle="Areg | Vegfa",x="X",y="Y"))

print(ggplot() +
        geom_sf(data=tissue_sf, fill="grey92", color="grey60", linewidth=0.4) +
        geom_point(data=tc_Areg,aes(x=x,y=y,color="Areg"),alpha=0.6,size=0.1) +
        geom_point(data=tc_Vegfa,aes(x=x,y=y,color="Vegfa"),alpha=0.5,size=0.1) +
        scale_color_manual(values=c("Areg"="#e74c3c","Vegfa"="#2980b9"),name=NULL) +
        #scale_size_continuous(range=c(0.3,2.5),name="log(UMI+1)") +
        coord_sf(datum=NA) + theme_demo() +
        labs(title="PNX7_M1_Sp",
             subtitle="Areg | Vegfa",x="X",y="Y"))+xlim(6000,6500)+ylim(11500,12000)

# ---- Section 2: Method comparison ----

cat("Section 2: Method comparison...\n")
pb <- list(D=1,
           diffusion_length = 10,
           lambda=0.3,
           production_rate=1,
           grid_resolution=256L,
           verbose=FALSE)

res_fd    <- do.call(estimate_concentration_field,
                     c(list(mask=tissue_mask,
                            transcript_coords=tc_Vegfa,
                            method="fd",
                            fd_solver="direct",
                            boundary_condition="dirichlet"),pb))
res_green <- do.call(estimate_concentration_field,
                     c(list(mask=tissue_mask,
                            transcript_coords=tc_Vegfa,
                            method="green"),pb))
res_kde   <- do.call(estimate_concentration_field,
                     c(list(mask=tissue_mask,
                            transcript_coords=tc_Vegfa,
                            method="kde",
                            kde_bandwidth=sqrt(1/0.3)),pb))
res_diff  <- res_fd; 
res_diff$field <- res_fd$field - res_green$field

print((plot_field(res_fd,tc_Vegfa,title="2a. fd (Dirichlet)",
                  subtitle="Exact BCs; recommended",
                  show_pts = F)+
         xlim(6000,6500)+
         ylim(11500,12000) +
         plot_field(res_green,tc_Vegfa,title="2b. green (FFT)",
                    subtitle="Infinite domain; no BCs",
                    show_pts = F)+
         xlim(6000,6500)+
         ylim(11500,12000)) /
        (plot_field(res_kde,tc_Vegfa,title="2c. kde",
                    subtitle="Bandwidth=L; phenomenological",
                    show_pts = F)+
           xlim(6000,6500)+
           ylim(11500,12000) +
           plot_field(res_diff,title="2d. fd \u2212 green",
                      subtitle="BCs suppress near-edge conc.",
                      symmetric=TRUE,
                      show_contours=FALSE,
                      show_pts=FALSE) +
        xlim(6000,6500)+
        ylim(11500,12000) +
        plot_annotation(title="Section 2: Method Comparison (Gene A, D=1, \u03bb=0.3)",
                        theme=theme(plot.title=element_text(face="bold",size=13)))))



