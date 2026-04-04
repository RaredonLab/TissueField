# ============================================================
# demo_fit_spatial_mask.R
#
# Basic demonstrations of fit_spatial_mask() across all
# three methods and their core tunable parameters.
#
# Run after sourcing: 
source("fit_spatial_mask.R")
#
# Required packages:
#   install.packages(c("sf","concaveman","MASS","isoband",
#                      "ggplot2","patchwork"))
#
# Sections:
#   1 — Circular cloud (sanity check)
#   2 — Elongated cloud (anisotropy)
#   3 — L-shaped cloud (non-convex; key diagnostic)
#   4 — Multi-cluster blobs
#   5 — Concavity parameter sweep
#   6 — Buffer and smoothing effects
#   7 — Validation: point containment table
# ============================================================

library(sf)
library(ggplot2)
library(patchwork)
source("fit_spatial_mask.R")
set.seed(42)

plot_mask <- function(mask, coords, title="", subtitle="",
                      fill="#4c9be8", border="#1a5fa8",
                      pt_col="#c0392b", pt_size=1.2, pt_alpha=0.55) {
  mc  <- as.data.frame(sf::st_coordinates(mask))
  grp <- if ("L3" %in% names(mc)) "L3" else if ("L2" %in% names(mc)) "L2" else "L1"
  ggplot() +
    geom_polygon(data=mc, aes(x=X,y=Y,group=.data[[grp]]),
                 fill=fill, alpha=0.2, color=border, linewidth=0.55) +
    geom_point(data=coords, aes(x=x,y=y), color=pt_col,
               size=pt_size, alpha=pt_alpha) +
    coord_equal() +
    labs(title=title, subtitle=subtitle, x="X", y="Y") +
    theme_minimal(base_size=9.5) +
    theme(plot.title=element_text(face="bold",size=9.5),
          plot.subtitle=element_text(size=8,color="grey40"),
          panel.grid=element_line(color="grey93"))
}

# ============================================================
# SECTION 1 — Circular cloud
# ============================================================
n <- 400
coords_circle <- data.frame(x=rnorm(n,0,1), y=rnorm(n,0,1))

m1a <- fit_spatial_mask(coords_circle, method="convex",  verbose=FALSE)
m1b <- fit_spatial_mask(coords_circle, method="concave", concavity=2, verbose=FALSE)
m1c <- fit_spatial_mask(coords_circle, method="kde", kde_threshold=0.02,
                         buffer_dist=0.1, verbose=FALSE)

print(
  plot_mask(m1a, coords_circle, title="1a. Circular — Convex",   subtitle="method='convex'") +
  plot_mask(m1b, coords_circle, title="1b. Circular — Concave",  subtitle="concavity=2") +
  plot_mask(m1c, coords_circle, title="1c. Circular — KDE",
            subtitle="threshold=0.02", fill="#27ae60", border="#1e8449") +
  plot_annotation(title="Section 1: Circular Cloud — Sanity Check",
    theme=theme(plot.title=element_text(face="bold",size=13)))
)

# ============================================================
# SECTION 2 — Elongated cloud
# ============================================================
coords_elong <- data.frame(x=rnorm(n,0,3), y=rnorm(n,0,0.6))

m2a <- fit_spatial_mask(coords_elong, method="convex",  verbose=FALSE)
m2b <- fit_spatial_mask(coords_elong, method="concave", verbose=FALSE)
m2c <- fit_spatial_mask(coords_elong, method="kde", kde_threshold=0.02,
                         buffer_dist=0.1, verbose=FALSE)

print(
  plot_mask(m2a, coords_elong, title="2a. Elongated — Convex") +
  plot_mask(m2b, coords_elong, title="2b. Elongated — Concave") +
  plot_mask(m2c, coords_elong, title="2c. Elongated — KDE",
            fill="#27ae60", border="#1e8449") +
  plot_annotation(title="Section 2: Elongated Cloud",
    theme=theme(plot.title=element_text(face="bold",size=13)))
)

# ============================================================
# SECTION 3 — L-shaped cloud (KEY test)
# ============================================================
coords_L <- rbind(
  data.frame(x=runif(200,0,5)+rnorm(200,0,0.05), y=runif(200,0,1)+rnorm(200,0,0.05)),
  data.frame(x=runif(200,0,1)+rnorm(200,0,0.05), y=runif(200,0,5)+rnorm(200,0,0.05)))

m3a <- fit_spatial_mask(coords_L, method="convex",  verbose=FALSE)
m3b <- fit_spatial_mask(coords_L, method="concave", concavity=1.5, verbose=FALSE)
m3c <- fit_spatial_mask(coords_L, method="concave", concavity=3.0, verbose=FALSE)

print(
  plot_mask(m3a, coords_L, title="3a. L-shape — Convex",
            subtitle="[includes empty corner]", fill="#e74c3c", border="#922b21") +
  plot_mask(m3b, coords_L, title="3b. Concave tight (1.5)",
            subtitle="[correctly excludes corner]") +
  plot_mask(m3c, coords_L, title="3c. Concave relaxed (3.0)",
            subtitle="[partially includes corner]", fill="#8e44ad", border="#6c3483") +
  plot_annotation(title="Section 3: L-Shaped Cloud — Non-Convex",
    subtitle="Lower concavity wraps more tightly; convex hull fails",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40")))
)

# ============================================================
# SECTION 4 — Multi-cluster
# ============================================================
mk_cluster <- function(n,cx,cy,sd=0.4) data.frame(x=rnorm(n,cx,sd),y=rnorm(n,cy,sd))
coords_multi <- rbind(mk_cluster(150,0,0), mk_cluster(150,4,1), mk_cluster(150,2,4))

m4a <- fit_spatial_mask(coords_multi, method="convex",  verbose=FALSE)
m4b <- fit_spatial_mask(coords_multi, method="concave", concavity=2, verbose=FALSE)
m4c <- fit_spatial_mask(coords_multi, method="kde", kde_threshold=0.01,
                         buffer_dist=0.15, verbose=FALSE)
m4d <- fit_spatial_mask(coords_multi, method="kde", kde_threshold=0.25,
                         buffer_dist=0.15, verbose=FALSE)

print(
  (plot_mask(m4a,coords_multi,title="4a. Multi — Convex") +
   plot_mask(m4b,coords_multi,title="4b. Multi — Concave")) /
  (plot_mask(m4c,coords_multi,title="4c. Multi — KDE inclusive",
             subtitle="threshold=0.01", fill="#27ae60", border="#1e8449") +
   plot_mask(m4d,coords_multi,title="4d. Multi — KDE tight",
             subtitle="threshold=0.25: separate masks", fill="#8e44ad", border="#6c3483")) +
  plot_annotation(title="Section 4: Multi-Cluster",
    theme=theme(plot.title=element_text(face="bold",size=13)))
)

# ============================================================
# SECTION 5 — Concavity sweep
# ============================================================
theta <- seq(0.2, pi-0.2, length.out=300)
coords_cres <- data.frame(
  x=cos(theta)*(1+rnorm(300,0,0.08))*3,
  y=sin(theta)*(1+rnorm(300,0,0.08))*3 - cos(theta)*1.5)

pal <- c("#e74c3c","#e67e22","#27ae60","#2980b9")
sw  <- mapply(function(cv,col) {
  m <- fit_spatial_mask(coords_cres, method="concave", concavity=cv, verbose=FALSE)
  plot_mask(m, coords_cres, title=paste0("concavity = ",cv), fill=col, border=col)
}, c(1.2,2,3,5), pal, SIMPLIFY=FALSE)

print(wrap_plots(sw, nrow=1) +
  plot_annotation(
    title="Section 5: Concavity Sweep — Crescent",
    subtitle="Lower = tighter | Higher = approaches convex",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40")))
)

# ============================================================
# SECTION 6 — Buffer and smoothing
# ============================================================
r   <- sqrt(runif(250,0.4^2,1^2)); ang <- runif(250,0,2*pi)
coords_ring <- data.frame(x=r*cos(ang)*5, y=r*sin(ang)*5)

m6a <- fit_spatial_mask(coords_ring, method="concave", concavity=2,
                         buffer_dist=0, smooth_mask=FALSE, verbose=FALSE)
m6b <- fit_spatial_mask(coords_ring, method="concave", concavity=2,
                         buffer_dist=0.3, smooth_mask=FALSE, verbose=FALSE)
m6c <- fit_spatial_mask(coords_ring, method="concave", concavity=2,
                         buffer_dist=0.3, smooth_mask=TRUE,
                         smooth_tolerance=0.15, verbose=FALSE)

print(
  plot_mask(m6a, coords_ring, title="6a. No buffer / no smooth") +
  plot_mask(m6b, coords_ring, title="6b. Buffer=0.3",
            fill="#27ae60", border="#1e8449") +
  plot_mask(m6c, coords_ring, title="6c. Buffer + smooth",
            fill="#8e44ad", border="#6c3483") +
  plot_annotation(title="Section 6: Buffer and Smoothing",
    theme=theme(plot.title=element_text(face="bold",size=13)))
)

# ============================================================
# SECTION 7 — Validation table
# ============================================================
cases <- list(
  list(label="Circle/convex",   mask=m1a, coords=coords_circle),
  list(label="Circle/concave",  mask=m1b, coords=coords_circle),
  list(label="Circle/kde",      mask=m1c, coords=coords_circle),
  list(label="Elongated/convex",mask=m2a, coords=coords_elong),
  list(label="Elongated/concave",mask=m2b,coords=coords_elong),
  list(label="Elongated/kde",   mask=m2c, coords=coords_elong),
  list(label="L-shape/convex",  mask=m3a, coords=coords_L),
  list(label="L-shape/conc1.5", mask=m3b, coords=coords_L),
  list(label="L-shape/conc3",   mask=m3c, coords=coords_L),
  list(label="Multi/convex",    mask=m4a, coords=coords_multi),
  list(label="Multi/concave",   mask=m4b, coords=coords_multi),
  list(label="Multi/kde-lo",    mask=m4c, coords=coords_multi),
  list(label="Multi/kde-hi",    mask=m4d, coords=coords_multi),
  list(label="Ring/raw",        mask=m6a, coords=coords_ring),
  list(label="Ring/buffer",     mask=m6b, coords=coords_ring),
  list(label="Ring/smooth",     mask=m6c, coords=coords_ring)
)

results <- do.call(rbind, lapply(cases, function(vc) {
  pts <- sf::st_geometry(sf::st_as_sf(vc$coords, coords=c("x","y")))
  ok  <- sf::st_within(pts, sf::st_union(vc$mask), sparse=FALSE)[,1L]
  n   <- nrow(vc$coords)
  data.frame(Case=vc$label, N=n, Enclosed=sum(ok),
             Pct=round(100*sum(ok)/n,2),
             Result=ifelse(sum(ok)==n,"\u2713 PASS","\u2717 FAIL"),
             stringsAsFactors=FALSE)
}))

cat("\n============================================================\n")
cat("  VALIDATION: Point Containment\n")
cat("============================================================\n")
print(results, row.names=FALSE)
cat("============================================================\n")
cat("demo_fit_spatial_mask.R complete.\n")

