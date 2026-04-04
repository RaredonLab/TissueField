# ============================================================
# demo_concentration_field.R
#
# Comprehensive demonstration of estimate_concentration_field().
#
# Run after sourcing:
#   source("fit_spatial_mask.R")
#   source("estimate_concentration_field.R")
#
# Required packages:
#   install.packages(c("sf","concaveman","MASS","isoband",
#                      "Matrix","ggplot2","patchwork"))
#
# Sections:
#   1  — Synthetic tissue + raw transcript layout
#   2  — Method comparison: fd / green / kde
#   3  — Diffusion length L sweep (most important parameter)
#   4  — Fixed L, varying D and lambda
#   5  — Dirichlet vs. Neumann boundary conditions
#   6  — UMI weighting vs. equal weighting
#   7  — External transcripts (include_external)
#   8  — Holey mask: vascular clearance topology
#   9  — Wide dynamic range + log1p visualisation
#   10 — Two-gene log2(A/B) ratio map
#   11 — 100k transcripts + solver timing comparison
#   12 — Quantitative field extraction
# ============================================================

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


# ============================================================
# SECTION 1 — Synthetic tissue
# ============================================================
cat("Section 1: Building tissue...\n")

in_kidney <- function(x,y,oa=9,ob=7,ncx=5,nr=3.8)
  (x/oa)^2+(y/ob)^2<1 & sqrt((x-ncx)^2+y^2)>=nr

cands      <- data.frame(x=runif(18000,-11,11),y=runif(18000,-9,9))
tissue_pts <- cands[in_kidney(cands$x,cands$y),][1:3000,]
tissue_mask <- fit_spatial_mask(tissue_pts, method="raster",
  raster_resolution=350L, raster_sigma=10, raster_threshold=0.15, verbose=FALSE)

mu_t <- sf::st_union(tissue_mask)

tc_A_cl <- data.frame(x=rnorm(350,-4.5,1.2), y=rnorm(350,3.8,0.9))
A_in    <- as.logical(sf::st_within(
  sf::st_as_sf(tc_A_cl,coords=c("x","y")),mu_t,sparse=FALSE)[,1L])
tc_A_cl <- tc_A_cl[A_in,]
tc_A_sc <- sample_in_mask(tissue_mask,50)
tc_A    <- rbind(tc_A_cl,tc_A_sc)
tc_A$umi <- c(rpois(nrow(tc_A_cl),15), rpois(nrow(tc_A_sc),2))

tc_B_all <- sample_in_mask(tissue_mask,2000)
prob_B   <- plogis(tc_B_all$x/3)
tc_B     <- tc_B_all[runif(nrow(tc_B_all))<prob_B,][1:500,]
tc_B$umi <- rpois(nrow(tc_B),5)

mdf <- as.data.frame(sf::st_coordinates(tissue_mask))
grp <- if ("L3" %in% names(mdf)) "L3" else if ("L2" %in% names(mdf)) "L2" else "L1"

print(ggplot() +
  geom_polygon(data=mdf,aes(x=X,y=Y,group=.data[[grp]]),
               fill="grey92",color="grey60",linewidth=0.4) +
  geom_point(data=tc_A,aes(x=x,y=y,color="Gene A",size=log1p(umi)),alpha=0.6) +
  geom_point(data=tc_B,aes(x=x,y=y,color="Gene B",size=log1p(umi)),alpha=0.5) +
  scale_color_manual(values=c("Gene A"="#e74c3c","Gene B"="#2980b9"),name=NULL) +
  scale_size_continuous(range=c(0.3,2.5),name="log(UMI+1)") +
  coord_equal() + theme_demo() +
  labs(title="Section 1: Synthetic kidney-bean tissue",
       subtitle="Gene A: focal upper-left | Gene B: broad right-biased",x="X",y="Y"))


# ============================================================
# SECTION 2 — Method comparison
# ============================================================
cat("Section 2: Method comparison...\n")
pb <- list(D=1,lambda=0.3,production_rate=1,grid_resolution=256L,verbose=FALSE)

res_fd    <- do.call(estimate_concentration_field,
  c(list(mask=tissue_mask,transcript_coords=tc_A,method="fd",
         fd_solver="direct",boundary_condition="dirichlet"),pb))
res_green <- do.call(estimate_concentration_field,
  c(list(mask=tissue_mask,transcript_coords=tc_A,method="green"),pb))
res_kde   <- do.call(estimate_concentration_field,
  c(list(mask=tissue_mask,transcript_coords=tc_A,method="kde",
         kde_bandwidth=sqrt(1/0.3)),pb))
res_diff  <- res_fd; res_diff$field <- res_fd$field - res_green$field

print((plot_field(res_fd,tc_A,title="2a. fd (Dirichlet)",subtitle="Exact BCs; recommended") +
       plot_field(res_green,tc_A,title="2b. green (FFT)",subtitle="Infinite domain; no BCs")) /
      (plot_field(res_kde,tc_A,title="2c. kde",subtitle="Bandwidth=L; phenomenological") +
       plot_field(res_diff,title="2d. fd \u2212 green",subtitle="BCs suppress near-edge conc.",
                  symmetric=TRUE,show_contours=FALSE,show_pts=FALSE)) +
  plot_annotation(title="Section 2: Method Comparison (Gene A, D=1, \u03bb=0.3)",
    theme=theme(plot.title=element_text(face="bold",size=13))))


# ============================================================
# SECTION 3 — Diffusion length sweep
# ============================================================
cat("Section 3: L sweep...\n")
L_vals <- c(0.5,1.5,4.0,10.0)
sw3 <- lapply(L_vals, function(Lv) {
  res <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
    D=1,diffusion_length=Lv,production_rate=1,method="fd",fd_solver="direct",
    grid_resolution=256L,verbose=FALSE)
  f <- res$field
  res$field <- (f-min(f,na.rm=TRUE))/diff(range(f,na.rm=TRUE))
  plot_field(res,tc_A,title=sprintf("L = %.1f",Lv),
             subtitle=sprintf("\u03bb \u2248 %.3f",1/Lv^2),
             fill_label="Norm. Conc.",pt_size=0.4) +
  scale_fill_viridis_c(option="magma",name="Norm.\nConc.",
                        limits=c(0,1),na.value="transparent")
})
print(wrap_plots(sw3,nrow=1) +
  plot_annotation(
    title="Section 3: Diffusion Length Sweep (normalised per panel)",
    subtitle="Small L \u2192 tight peaks | Large L \u2192 tissue-wide field",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40"))))

peaks <- sapply(L_vals, function(Lv) {
  res <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
    D=1,diffusion_length=Lv,method="fd",fd_solver="direct",
    grid_resolution=256L,verbose=FALSE)
  max(res$field,na.rm=TRUE)
})
cat("\n  Peak concentration vs. L:\n")
print(data.frame(L=L_vals, peak=round(peaks,4)))
cat("  Start with L ~ 10-20% of inter-cluster distance.\n")


# ============================================================
# SECTION 4 — Fixed L=2, varying D and lambda
# ============================================================
cat("\nSection 4: Fixed L=2, varying D...\n")
DL <- list(list(D=0.1,lam=0.025,label="D=0.1, \u03bb=0.025"),
           list(D=1,  lam=0.25, label="D=1,   \u03bb=0.25 "),
           list(D=5,  lam=1.25, label="D=5,   \u03bb=1.25 "),
           list(D=20, lam=5,    label="D=20,  \u03bb=5    "))
dl_pl <- lapply(DL, function(co) {
  res <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
    D=co$D,lambda=co$lam,production_rate=1,method="fd",fd_solver="direct",
    grid_resolution=256L,verbose=FALSE)
  pk <- round(max(res$field,na.rm=TRUE),5)
  plot_field(res,tc_A,title=co$label,
             subtitle=sprintf("L=2 fixed | peak=%.5f",pk),pt_size=0.3)
})
print(wrap_plots(dl_pl,nrow=1) +
  plot_annotation(
    title="Section 4: Fixed L=2, Varying D and \u03bb",
    subtitle="Shape IDENTICAL (same L). Amplitude \u221d 1/D.",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40"))))


# ============================================================
# SECTION 5 — Boundary conditions
# ============================================================
cat("\nSection 5: Boundary conditions...\n")
res_dir <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
  D=1,lambda=0.2,method="fd",fd_solver="direct",boundary_condition="dirichlet",
  grid_resolution=256L,verbose=FALSE)
res_neu <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
  D=1,lambda=0.2,method="fd",fd_solver="direct",boundary_condition="neumann",
  grid_resolution=256L,verbose=FALSE)
res_bcd <- res_dir; res_bcd$field <- res_neu$field - res_dir$field

sl <- function(res,lbl) {
  df <- field_to_df(res)
  s  <- df[abs(df$y)<res$hy*1.5&!is.na(df$field),]
  s  <- s[order(s$x),]; s$BC <- lbl; s
}
sl_all <- rbind(sl(res_dir,"Dirichlet"),sl(res_neu,"Neumann"))

p5d <- ggplot(sl_all,aes(x=x,y=field,color=BC)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values=c("Dirichlet"="#e74c3c","Neumann"="#2980b9")) +
  labs(title="5d. Cross-section y\u22480",
       subtitle="Dirichlet drops to 0; Neumann stays elevated",
       x="X",y="Conc.",color="BC") + theme_demo()

print((plot_field(res_dir,tc_A,title="5a. Dirichlet (C=0 at boundary)",
                  subtitle="D=1, \u03bb=0.2, L\u22482.24") +
       plot_field(res_neu,tc_A,title="5b. Neumann (dC/dn=0)",
                  subtitle="Molecules reflected",palette="plasma")) /
      (plot_field(res_bcd,title="5c. Neumann \u2212 Dirichlet",
                  subtitle="Neumann retains more near-boundary",
                  symmetric=TRUE,show_contours=FALSE,show_pts=FALSE) + p5d) +
  plot_annotation(title="Section 5: Boundary Condition Comparison",
    theme=theme(plot.title=element_text(face="bold",size=13))))


# ============================================================
# SECTION 6 — UMI weighting
# ============================================================
cat("\nSection 6: UMI weighting...\n")
res_uw <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
  D=1,lambda=0.3,method="fd",fd_solver="direct",grid_resolution=256L,
  weight_col=NULL,verbose=FALSE)
res_wt <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_A,
  D=1,lambda=0.3,method="fd",fd_solver="direct",grid_resolution=256L,
  weight_col="umi",verbose=FALSE)
nrm <- function(r) { f<-r$field; r$field<-(f-min(f,na.rm=TRUE))/diff(range(f,na.rm=TRUE)); r }
print(
  plot_field(nrm(res_uw),tc_A,title="6a. Unweighted",
             subtitle="All transcripts equal",fill_label="Norm. Conc.") +
  plot_field(nrm(res_wt),tc_A,title="6b. UMI-weighted",
             subtitle="High-UMI cluster amplified",palette="plasma",
             fill_label="Norm. Conc.") +
  plot_annotation(title="Section 6: UMI Weighting",
    theme=theme(plot.title=element_text(face="bold",size=13))))


# ============================================================
# SECTION 7 — External transcripts
# ============================================================
cat("\nSection 7: External transcripts...\n")
tc_ext <- data.frame(x=rnorm(80,6.5,0.4),y=rnorm(80,0,0.5),umi=rpois(80,10))
ext_in <- as.logical(sf::st_within(
  sf::st_as_sf(tc_ext,coords=c("x","y")),mu_t,sparse=FALSE)[,1L])
tc_ext <- tc_ext[!ext_in,]
tc_Ax  <- rbind(tc_A,tc_ext)
res_ne <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_Ax,
  D=1,lambda=0.15,method="fd",fd_solver="direct",grid_resolution=256L,
  include_external=FALSE,verbose=FALSE)
res_we <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_Ax,
  D=1,lambda=0.15,method="fd",fd_solver="direct",grid_resolution=256L,
  include_external=TRUE,verbose=FALSE)
ext_pt_layer <- geom_point(data=tc_ext,aes(x=x,y=y),color="orange",
                            size=0.9,alpha=0.8,inherit.aes=FALSE)
print(
  (plot_field(res_ne,tc_A,title="7a. include_external=FALSE") + ext_pt_layer) +
  (plot_field(res_we,tc_A,title="7b. include_external=TRUE",palette="plasma") + ext_pt_layer) +
  plot_annotation(title="Section 7: External Transcripts (orange=external)",
    theme=theme(plot.title=element_text(face="bold",size=13))))


# ============================================================
# SECTION 8 — Holey mask (vasculature)
# ============================================================
cat("\nSection 8: Holey mask...\n")
vd <- list(list(cx=-2,cy=1.5,r=1.2),list(cx=2.5,cy=-2,r=1),list(cx=-4.5,cy=-2.5,r=0.8))
in_v <- function(x,y) Reduce(`|`,lapply(vd,function(v) sqrt((x-v$cx)^2+(y-v$cy)^2)<v$r))

n_h  <- 3000; cnd_h <- data.frame(x=runif(n_h*8,-11,11),y=runif(n_h*8,-9,9))
kp_h <- in_kidney(cnd_h$x,cnd_h$y) & !in_v(cnd_h$x,cnd_h$y)
pts_h <- cnd_h[kp_h,][1:n_h,]
holey_mask <- fit_spatial_mask(pts_h,method="raster",
  raster_resolution=400L,raster_sigma=7,raster_threshold=0.2,verbose=FALSE)
tc_h <- sample_in_mask(holey_mask,400)

res_hfd  <- estimate_concentration_field(mask=holey_mask,transcript_coords=tc_h,
  D=1,lambda=0.15,method="fd",fd_solver="direct",grid_resolution=300L,verbose=FALSE)
res_htight <- estimate_concentration_field(mask=holey_mask,transcript_coords=tc_h,
  D=1,lambda=1.0,method="fd",fd_solver="direct",grid_resolution=300L,verbose=FALSE)
res_hgrn <- estimate_concentration_field(mask=holey_mask,transcript_coords=tc_h,
  D=1,lambda=0.15,method="green",grid_resolution=300L,verbose=FALSE)
res_hdiff <- res_hfd; res_hdiff$field <- res_hfd$field - res_hgrn$field

print((plot_field(res_hfd,tc_h,title="8a. fd | L\u22482.58",
                  subtitle="Vessels = clearance sinks",pt_size=0.3) +
       plot_field(res_htight,tc_h,title="8b. fd | L=1 (tight)",
                  subtitle="Steeper depletion near vessels",palette="inferno",pt_size=0.3)) /
      (plot_field(res_hgrn,tc_h,title="8c. green | holes invisible",
                  subtitle="FFT ignores topology",palette="plasma",pt_size=0.3) +
       plot_field(res_hdiff,title="8d. fd \u2212 green",
                  subtitle="fd correctly zeros lumens; green does not",
                  symmetric=TRUE,show_contours=FALSE,show_pts=FALSE)) +
  plot_annotation(title="Section 8: Holey Mask — Vascular Clearance",
    subtitle="fd zeros vessel lumens; green fills them — use fd for holey domains",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40"))))


# ============================================================
# SECTION 9 — Wide dynamic range + log scale
# ============================================================
cat("\nSection 9: Wide dynamic range...\n")
tc_w <- rbind(data.frame(x=rnorm(20,-5,0.3),y=rnorm(20,4,0.3),umi=rpois(20,200)),
              data.frame(x=rnorm(300,0,4),   y=rnorm(300,0,3), umi=rpois(300,1)))
tc_w_in <- as.logical(sf::st_within(
  sf::st_as_sf(tc_w[,c("x","y")],coords=c("x","y")),mu_t,sparse=FALSE)[,1L])
tc_w <- tc_w[tc_w_in,]
res_w <- estimate_concentration_field(mask=tissue_mask,transcript_coords=tc_w,
  D=1,lambda=0.5,method="fd",fd_solver="direct",grid_resolution=256L,
  weight_col="umi",verbose=FALSE)
print(
  plot_field(res_w,tc_w,title="9a. Linear",
             subtitle="Only cluster visible",pt_color="white") +
  plot_field(res_w,tc_w,log_scale=TRUE,title="9b. log\u2081\u208a",
             subtitle="Background + cluster both visible",
             palette="viridis",pt_color="white") +
  plot_annotation(title="Section 9: Wide Dynamic Range",
    subtitle="Use log_scale=TRUE when max/min UMI > 20",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40"))))


# ============================================================
# SECTION 10 — Two genes + ratio map
# ============================================================
cat("\nSection 10: Two-gene ratio...\n")
ph <- list(D=1,lambda=0.25,production_rate=1,method="fd",fd_solver="direct",
           grid_resolution=256L,verbose=FALSE)
res_gA <- do.call(estimate_concentration_field,c(list(mask=tissue_mask,transcript_coords=tc_A),ph))
res_gB <- do.call(estimate_concentration_field,c(list(mask=tissue_mask,transcript_coords=tc_B),ph))
res_rt <- res_gA; res_rt$field <- log2((res_gA$field+1e-6)/(res_gB$field+1e-6))

print((plot_field(res_gA,tc_A,title="10a. Gene A") +
       plot_field(res_gB,tc_B,title="10b. Gene B",palette="viridis")) /
      (plot_field(res_rt,title="10c. log\u2082(A/B) ratio",
                  subtitle="Yellow=A territory | Blue=B | White=balanced",
                  symmetric=TRUE,show_contours=TRUE,n_contours=5,
                  show_pts=FALSE,fill_label="log\u2082(A/B)") +
       geom_point(data=tc_A,aes(x=x,y=y),color="#e74c3c",size=0.3,alpha=0.4,inherit.aes=FALSE) +
       geom_point(data=tc_B,aes(x=x,y=y),color="#27ae60",size=0.3,alpha=0.4,inherit.aes=FALSE) |
       ggplot(field_to_df(res_rt),aes(x=field,fill=ifelse(x<0,"Left","Right"))) +
       geom_density(alpha=0.55) +
       geom_vline(xintercept=0,linetype="dashed",color="grey40") +
       scale_fill_manual(values=c("Left"="#c0392b","Right"="#2980b9"),name="Half") +
       labs(title="10d. log\u2082(A/B) density",x="log\u2082(A/B)",y="Density") + theme_demo()) +
  plot_annotation(title="Section 10: Two-Gene Ratio Map",
    theme=theme(plot.title=element_text(face="bold",size=13))))


# ============================================================
# SECTION 11 — 100k transcripts + solver timing
# ============================================================
cat("\nSection 11: 100k scale test...\n")
tc_100k <- sample_in_mask(tissue_mask, 100000)
cat(sprintf("  Sampled %d transcripts.\n", nrow(tc_100k)))

cat("  fd/direct N=256...\n")
t1 <- system.time(r1 <- estimate_concentration_field(
  mask=tissue_mask,transcript_coords=tc_100k,D=1,lambda=0.3,
  method="fd",fd_solver="direct",grid_resolution=256L,verbose=TRUE))
cat("  fd/SOR N=512...\n")
t2 <- system.time(r2 <- estimate_concentration_field(
  mask=tissue_mask,transcript_coords=tc_100k,D=1,lambda=0.3,
  method="fd",fd_solver="iterative",sor_omega=1.75,grid_resolution=512L,verbose=TRUE))
cat("  green/FFT N=512...\n")
t3 <- system.time(r3 <- estimate_concentration_field(
  mask=tissue_mask,transcript_coords=tc_100k,D=1,lambda=0.3,
  method="green",grid_resolution=512L,verbose=TRUE))

cat("\n")
print(data.frame(
  Solver=c("fd/direct (N=256)","fd/SOR (N=512)","green/FFT (N=512)"),
  Time_s=round(c(t1["elapsed"],t2["elapsed"],t3["elapsed"]),2),
  Has_BCs=c("Yes","Yes","No"), Handles_holes=c("Yes","Yes","No")),
  row.names=FALSE)

sub <- tc_100k[sample(nrow(tc_100k),2000),]
print(
  plot_field(r1,sub,title=sprintf("11a. fd/direct (%.1fs)",t1["elapsed"]),
             subtitle="N=256 | sparse LU",pt_size=0.15,pt_alpha=0.2) +
  plot_field(r2,sub,title=sprintf("11b. fd/SOR (%.1fs)",t2["elapsed"]),
             subtitle="N=512 | iterative",palette="plasma",pt_size=0.15,pt_alpha=0.2) +
  plot_field(r3,sub,title=sprintf("11c. green (%.1fs)",t3["elapsed"]),
             subtitle="N=512 | FFT",palette="viridis",pt_size=0.15,pt_alpha=0.2) +
  plot_annotation(title="Section 11: 100k Transcripts — Solver Comparison",
    subtitle="2,000 of 100,000 points shown",
    theme=theme(plot.title=element_text(face="bold",size=13),
                plot.subtitle=element_text(size=9,color="grey40"))))


# ============================================================
# SECTION 12 — Quantitative extraction
# ============================================================
cat("\nSection 12: Quantitative extraction...\n")
res_q <- res_fd

interpolate_field <- function(result, query_pts) {
  xg<-result$x; yg<-result$y; f<-result$field
  vapply(seq_len(nrow(query_pts)), function(i) {
    qx<-query_pts$x[i]; qy<-query_pts$y[i]
    xi<-findInterval(qx,xg); yi<-findInterval(qy,yg)
    if (xi<1||xi>=length(xg)||yi<1||yi>=length(yg)) return(NA_real_)
    wx<-(qx-xg[xi])/(xg[xi+1]-xg[xi]); wy<-(qy-yg[yi])/(yg[yi+1]-yg[yi])
    (1-wx)*(1-wy)*f[xi,yi]+wx*(1-wy)*f[xi+1,yi]+(1-wx)*wy*f[xi,yi+1]+wx*wy*f[xi+1,yi+1]
  }, numeric(1))
}

queries <- data.frame(
  label=c("Cluster core","Cluster edge","Tissue centre","Far right","Near notch"),
  x=c(-4.5,-3.0,0.0,5.0,3.0), y=c(3.8,2.5,0.0,0.0,0.0))
queries$concentration <- interpolate_field(res_q, queries)

cat("\n  12a. Concentration at query points:\n")
print(queries[,c("label","x","y","concentration")], row.names=FALSE)

df_q <- field_to_df(res_q)
cat(sprintf("\n  12b. Mean conc. by half:\n    Left : %.5f\n    Right: %.5f\n",
    mean(df_q$field[df_q$x< 0 & !is.na(df_q$field)]),
    mean(df_q$field[df_q$x>=0 & !is.na(df_q$field)])))

pk  <- which.max(res_q$field); N <- length(res_q$x)
cat(sprintf("\n  12c. Peak: x=%.2f, y=%.2f, value=%.6f\n",
    res_q$x[((pk-1L)%%N)+1L], res_q$y[((pk-1L)%/%N)+1L], res_q$field[pk]))

total   <- sum(res_q$field[res_q$mask],na.rm=TRUE)*res_q$hx*res_q$hy
theory  <- nrow(tc_A)*1.0/0.3
cat(sprintf("\n  12d. Field integral   = %.4f\n", total))
cat(sprintf("       Theory (n*r/lam) = %.4f\n", theory))
cat(sprintf("       Agreement        = %.1f%%\n", 100*(1-abs(total-theory)/theory)))

prof <- df_q[abs(df_q$y)<res_q$hy*1.5 & !is.na(df_q$field),]
print(ggplot(prof[order(prof$x),], aes(x=x,y=field)) +
  geom_line(color="#e74c3c",linewidth=0.8) +
  geom_area(fill="#e74c3c",alpha=0.15) +
  geom_vline(xintercept=0,linetype="dashed",color="grey50") +
  labs(title="12e. Concentration profile along y\u22480",x="X",y="Conc.") + theme_demo())

cat("\n",strrep("=",60),"\n")
cat("  demo_concentration_field.R complete.\n\n")
cat("  Parameter guide:\n")
cat("  Method  : 'fd' default. 'green' fast interior. 'kde' baseline.\n")
cat("  L sweep : 0.5, 1, 2, 4, 8 (10-20% of inter-cluster dist).\n")
cat("  BC      : 'dirichlet' tissue sections. 'neumann' sealed.\n")
cat("  N       : 256 explore. 512 final. SOR for N > 400.\n")
cat(strrep("=",60),"\n")
