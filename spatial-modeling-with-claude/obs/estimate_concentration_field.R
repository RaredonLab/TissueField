# ============================================================
# estimate_concentration_field.R
#
# Models steady-state continuous molecular concentration over a
# spatial mask from individual mRNA transcript point sources.
#
# Physical model:  D * nabla^2 C(x) - lambda * C(x) + s(x) = 0
#
# Three solvers:
#   "fd"    — finite-difference sparse linear system (recommended)
#   "green" — FFT convolution with 2D Green's function K0(kappa*r)
#   "kde"   — Gaussian kernel smoothing (phenomenological baseline)
#
# install.packages(c("sf","Matrix","ggplot2"))
# ============================================================

.gauss_kernel_1d <- function(sigma) {
  r <- max(1L,ceiling(3*sigma)); k <- exp(-(seq(-r,r))^2/(2*sigma^2)); k/sum(k)
}

.sep_gauss_smooth_2d <- function(mat, sigma_x, sigma_y) {
  kx <- .gauss_kernel_1d(sigma_x); ky <- .gauss_kernel_1d(sigma_y)
  out <- apply(mat,2L,function(col){s<-stats::filter(col,kx,method="convolution",sides=2L);s[is.na(s)]<-0;as.numeric(s)})
  t(apply(out,1L,function(row){s<-stats::filter(row,ky,method="convolution",sides=2L);s[is.na(s)]<-0;as.numeric(s)}))
}

.rasterize_mask <- function(gc_sf, mask_union, N, n_cores) {
  n_cells <- N*N
  if (n_cores>1L && .Platform$OS.type=="unix") {
    chunk_size <- ceiling(n_cells/n_cores)
    chunks     <- split(seq_len(n_cells),ceiling(seq_len(n_cells)/chunk_size))
    res_list   <- parallel::mclapply(chunks,function(idx)
      as.logical(sf::st_within(gc_sf[idx,],mask_union,sparse=FALSE)[,1L]),mc.cores=n_cores)
    in_msk <- unlist(res_list)
  } else {
    in_msk <- as.logical(sf::st_within(gc_sf,mask_union,sparse=FALSE)[,1L])
  }
  matrix(in_msk,nrow=N,ncol=N)
}

estimate_concentration_field <- function(
  mask,
  transcript_coords,
  D                  = 1.0,
  lambda             = 0.1,
  production_rate    = 1.0,
  diffusion_length   = NULL,
  weight_col         = NULL,
  method             = "fd",
  grid_resolution    = 256L,
  boundary_condition = "dirichlet",
  fd_solver          = "direct",
  sor_omega          = 1.7,
  sor_tol            = 1e-6,
  sor_max_iter       = 1500L,
  green_r_min_factor = 0.5,
  kde_bandwidth      = NULL,
  include_external   = FALSE,
  normalize          = FALSE,
  log_transform      = FALSE,
  clip_negative      = TRUE,
  n_cores            = 1L,
  return_sources     = TRUE,
  plot               = FALSE,
  verbose            = TRUE
) {
  t0 <- proc.time()["elapsed"]

  # 0. Packages
  req  <- "sf"; if (method=="fd"&&fd_solver=="direct") req <- c(req,"Matrix")
  miss <- req[!sapply(req,requireNamespace,quietly=TRUE)]
  if (length(miss)>0L) stop("Install: install.packages(c(",paste0('"',miss,'"',collapse=","),"))")

  # 1. Physics
  if (!is.null(diffusion_length)) {
    L <- as.numeric(diffusion_length); kappa2 <- 1/L^2; kappa <- 1/L
  } else {
    if (lambda<=0 && method!="kde")
      stop("`lambda` must be > 0 for a finite steady-state field.")
    if (lambda>0) { L <- sqrt(D/lambda); kappa2 <- lambda/D; kappa <- sqrt(kappa2)
    } else { L <- Inf; kappa2 <- 0; kappa <- 0 }
  }
  N <- as.integer(grid_resolution); nx_c <- as.integer(n_cores)

  if (verbose) {
    cat("============================================================\n")
    cat("  estimate_concentration_field\n")
    cat("============================================================\n")
    cat(sprintf("  Method: %s | D: %g | lambda: %g | L: %g | N: %d x %d\n",
                method,D,lambda,L,N,N))
    if (method=="fd") cat(sprintf("  BC: %s | solver: %s\n",boundary_condition,fd_solver))
    cat("------------------------------------------------------------\n")
  }

  # 2. Parse transcripts
  if (is.matrix(transcript_coords)) transcript_coords <- as.data.frame(transcript_coords)
  if (!all(c("x","y") %in% names(transcript_coords))) names(transcript_coords)[1:2] <- c("x","y")
  tc      <- transcript_coords[complete.cases(transcript_coords[,c("x","y")]),]
  weights <- if (!is.null(weight_col)) tc[[weight_col]]*production_rate else
    rep(production_rate,nrow(tc))

  # 3. Filter by mask
  mu_sf  <- sf::st_union(mask)
  tc_sf  <- sf::st_as_sf(tc[,c("x","y")],coords=c("x","y"),crs=sf::st_crs(mask))
  in_msk <- as.logical(sf::st_within(tc_sf,mu_sf,sparse=FALSE)[,1L])
  n_tot  <- nrow(tc); n_ext <- sum(!in_msk)
  if (!include_external) { tc <- tc[in_msk,]; weights <- weights[in_msk] }
  n_src <- nrow(tc)
  if (n_src==0L) stop("No transcripts available. Check coords or set include_external=TRUE.")
  if (verbose) cat(sprintf("  Transcripts: %d used (%d external %s)\n",
                            n_src, n_ext, if(include_external)"included" else "excluded"))

  # 4. Grid
  bb   <- sf::st_bbox(mu_sf); pad <- 0.03
  xpad <- (bb["xmax"]-bb["xmin"])*pad; ypad <- (bb["ymax"]-bb["ymin"])*pad
  x_breaks  <- seq(bb["xmin"]-xpad,bb["xmax"]+xpad,length.out=N+1L)
  y_breaks  <- seq(bb["ymin"]-ypad,bb["ymax"]+ypad,length.out=N+1L)
  x_centers <- (x_breaks[-1L]+x_breaks[-(N+1L)])*0.5
  y_centers <- (y_breaks[-1L]+y_breaks[-(N+1L)])*0.5
  hx <- x_centers[2L]-x_centers[1L]; hy <- y_centers[2L]-y_centers[1L]; h <- (hx+hy)*0.5

  # 5. Rasterize mask
  if (verbose) cat(sprintf("  Rasterizing %dx%d grid...\n",N,N))
  gdf     <- expand.grid(x=x_centers,y=y_centers)
  gc_sf   <- sf::st_as_sf(gdf,coords=c("x","y"),crs=sf::st_crs(mask))
  mask_mat <- .rasterize_mask(gc_sf,mu_sf,N,nx_c)
  n_inside <- sum(mask_mat)
  if (verbose) cat(sprintf("  Inside-mask cells: %d (%.1f%%)\n",n_inside,100*n_inside/N^2))

  # 6. Bin sources
  xi_s <- pmax(1L,pmin(N,findInterval(tc$x,x_breaks,rightmost.closed=TRUE)))
  yi_s <- pmax(1L,pmin(N,findInterval(tc$y,y_breaks,rightmost.closed=TRUE)))
  lin_s <- xi_s+(yi_s-1L)*N
  src_vec <- numeric(N*N)
  wt_agg  <- rowsum(matrix(weights,ncol=1L),lin_s)
  src_vec[as.integer(rownames(wt_agg))] <- wt_agg[,1L]
  source_mat <- matrix(src_vec,nrow=N,ncol=N)

  # 7. Solve
  t_solve <- proc.time()["elapsed"]

  field_mat <- switch(method,

    "fd" = {
      inside_lin <- which(mask_mat); K <- length(inside_lin)
      cell2sys <- integer(N*N); cell2sys[inside_lin] <- seq_len(K)
      xi_in <- ((inside_lin-1L)%%N)+1L; yi_in <- ((inside_lin-1L)%/%N)+1L
      nbr_xm <- ifelse(xi_in>1L,inside_lin-1L,0L)
      nbr_xp <- ifelse(xi_in<N, inside_lin+1L,0L)
      nbr_ym <- ifelse(yi_in>1L,inside_lin-N, 0L)
      nbr_yp <- ifelse(yi_in<N, inside_lin+N, 0L)
      nx_w <- 1/hx^2; ny_w <- 1/hy^2
      is_in_mask <- function(nb){ok<-nb>0L&nb<=N*N;res<-logical(length(nb));res[ok]<-mask_mat[nb[ok]];res}
      nbr_list <- list(nbr_xm,nbr_xp,nbr_ym,nbr_yp); wt_list <- c(nx_w,nx_w,ny_w,ny_w)

      if (fd_solver=="direct") {
        diag_vals <- rep(-(2*nx_w+2*ny_w+kappa2),K)
        off_rows <- integer(0); off_cols <- integer(0); off_vals <- numeric(0)
        for (d in seq_along(nbr_list)) {
          nb <- nbr_list[[d]]; wt <- wt_list[[d]]
          is_in <- is_in_mask(nb); is_out <- (nb>0L)&!is_in; is_oob <- (nb==0L)
          sel <- which(is_in)
          if (length(sel)>0L){off_rows<-c(off_rows,sel);off_cols<-c(off_cols,cell2sys[nb[sel]]);off_vals<-c(off_vals,rep(wt,length(sel)))}
          if (boundary_condition=="neumann") {
            nc2 <- which(is_out|is_oob); diag_vals[nc2] <- diag_vals[nc2]+wt
          }
        }
        A <- Matrix::sparseMatrix(i=c(seq_len(K),off_rows),j=c(seq_len(K),off_cols),
                                   x=c(diag_vals,off_vals),dims=c(K,K))
        b_vec <- -source_mat[inside_lin]/(D*hx*hy)
        if (verbose) cat(sprintf("  Sparse system: %dx%d, %d nnz\n",K,K,length(A@x)))
        c_vec <- as.numeric(Matrix::solve(A,b_vec))
        out <- matrix(NA_real_,N,N); out[inside_lin] <- c_vec; out

      } else {
        if (boundary_condition=="neumann") message("SOR uses Dirichlet BC only.")
        xi_grid <- matrix(rep(seq_len(N),N),N,N)
        yi_grid <- matrix(rep(seq_len(N),each=N),N,N)
        is_red  <- ((xi_grid+yi_grid)%%2L)==0L
        ins_red <- mask_mat& is_red; ins_blk <- mask_mat&!is_red
        C <- matrix(0,N,N); denom <- 2*nx_w+2*ny_w+kappa2
        src_dens <- source_mat/(D*hx*hy); converged <- FALSE
        for (iter in seq_len(sor_max_iter)) {
          C_prev <- C
          Cxm <- rbind(matrix(0,1,N),C[1:(N-1),]); Cxp <- rbind(C[2:N,],matrix(0,1,N))
          Cym <- cbind(matrix(0,N,1),C[,1:(N-1)]); Cyp <- cbind(C[,2:N],matrix(0,N,1))
          Cgs <- (nx_w*(Cxm+Cxp)+ny_w*(Cym+Cyp)+src_dens)/denom
          C[ins_red] <- (1-sor_omega)*C[ins_red]+sor_omega*Cgs[ins_red]; C[!mask_mat]<-0
          Cxm <- rbind(matrix(0,1,N),C[1:(N-1),]); Cxp <- rbind(C[2:N,],matrix(0,1,N))
          Cym <- cbind(matrix(0,N,1),C[,1:(N-1)]); Cyp <- cbind(C[,2:N],matrix(0,N,1))
          Cgs <- (nx_w*(Cxm+Cxp)+ny_w*(Cym+Cyp)+src_dens)/denom
          C[ins_blk] <- (1-sor_omega)*C[ins_blk]+sor_omega*Cgs[ins_blk]; C[!mask_mat]<-0
          delta <- max(abs(C[mask_mat]-C_prev[mask_mat]))
          if (is.finite(delta)&&delta<sor_tol){if(verbose)cat(sprintf("  SOR converged: iter=%d, delta=%.2e\n",iter,delta));converged<-TRUE;break}
        }
        if (!converged) warning("SOR did not converge in ",sor_max_iter," iterations.")
        out <- matrix(NA_real_,N,N); out[mask_mat] <- C[mask_mat]; out
      }
    },

    "green" = {
      if (kappa<=0) stop("method='green' requires lambda>0.")
      Np   <- 2L*N
      xi_k <- c(0L:(N-1L),(-N):(-1L)); yi_k <- c(0L:(N-1L),(-N):(-1L))
      dx_mat <- outer(xi_k,rep(1L,Np))*hx; dy_mat <- outer(rep(1L,Np),yi_k)*hy
      r_mat  <- sqrt(dx_mat^2+dy_mat^2); r_min <- green_r_min_factor*h
      r_mat[r_mat<r_min] <- r_min
      G_mat  <- besselK(kappa*r_mat,nu=0)*(hx*hy)/(2*pi*D)
      G_fft  <- fft(G_mat)
      S_pad  <- matrix(0,Np,Np); S_pad[1L:N,1L:N] <- source_mat/(hx*hy)
      C_conv <- Re(fft(G_fft*fft(S_pad),inverse=TRUE))/Np^2
      out <- matrix(NA_real_,N,N); out[mask_mat] <- C_conv[1L:N,1L:N][mask_mat]; out
    },

    "kde" = {
      bw <- if (!is.null(kde_bandwidth)) as.numeric(kde_bandwidth) else
        if (is.finite(L)) L else min(diff(range(x_centers)),diff(range(y_centers)))*0.05
      if (verbose) cat(sprintf("  KDE bandwidth: %g\n",bw))
      sm  <- .sep_gauss_smooth_2d(source_mat,bw/hx,bw/hy)
      out <- matrix(NA_real_,N,N); out[mask_mat] <- sm[mask_mat]; out
    },

    stop("Unknown method '",method,"'. Choose: fd, green, kde.")
  )

  t_solve_end <- proc.time()["elapsed"]

  # 8. Post-processing
  field_mat[!mask_mat] <- NA_real_
  if (clip_negative) field_mat[!is.na(field_mat)&field_mat<0] <- 0
  if (log_transform)  field_mat <- log1p(field_mat)
  if (normalize) {
    fmn <- min(field_mat,na.rm=TRUE); fmx <- max(field_mat,na.rm=TRUE)
    if (fmx>fmn) field_mat <- (field_mat-fmn)/(fmx-fmn)
  }

  t_total <- proc.time()["elapsed"]-t0
  if (verbose) {
    cat(sprintf("  Solve: %.2fs | Total: %.2fs | Range: [%.4g, %.4g]\n",
                t_solve_end-t_solve, t_total,
                min(field_mat,na.rm=TRUE),max(field_mat,na.rm=TRUE)))
    cat("============================================================\n")
  }

  # 9. Plot
  if (plot && requireNamespace("ggplot2",quietly=TRUE)) {
    gdf <- expand.grid(x=x_centers,y=y_centers); gdf$field <- as.vector(field_mat)
    gdf <- gdf[!is.na(gdf$field),]
    print(ggplot2::ggplot(gdf,ggplot2::aes(x=x,y=y,fill=field))+
      ggplot2::geom_raster(interpolate=TRUE)+
      ggplot2::scale_fill_viridis_c(option="magma",name="Conc.",na.value="transparent")+
      ggplot2::coord_equal()+ggplot2::theme_minimal(base_size=11)+
      ggplot2::labs(title=paste0("Concentration | method='",method,"'"),x="X",y="Y"))
  }

  # 10. Return
  list(
    field   = field_mat,
    x       = x_centers, y = y_centers, hx = hx, hy = hy,
    mask    = mask_mat,
    sources = if (return_sources) source_mat else NULL,
    params  = list(D=D,lambda=lambda,kappa2=kappa2,diffusion_length=L,
                   production_rate=production_rate,method=method,
                   grid_resolution=N,boundary_condition=if(method=="fd")boundary_condition else NA,
                   fd_solver=if(method=="fd")fd_solver else NA,
                   include_external=include_external,log_transform=log_transform,normalize=normalize),
    diagnostics = list(n_transcripts_total=n_tot,n_transcripts_used=n_src,
                       n_transcripts_external=n_ext,n_cells_inside_mask=n_inside,
                       elapsed_total_s=t_total,elapsed_solve_s=t_solve_end-t_solve)
  )
}

# ============================================================
# UTILITY FUNCTIONS
# ============================================================

field_to_df <- function(result, inside_only=TRUE) {
  df <- expand.grid(x=result$x,y=result$y)
  df$field  <- as.vector(result$field)
  df$inside <- as.vector(result$mask)
  if (!is.null(result$sources)) df$source <- as.vector(result$sources)
  if (inside_only) df <- df[df$inside,]
  df
}

plot_concentration_field <- function(
    result, transcript_coords=NULL,
    show_sources=TRUE, show_contours=FALSE, n_contours=8L,
    palette="magma", interpolate=TRUE, log_scale=FALSE,
    pt_size=0.4, pt_alpha=0.4, pt_color="#00e5ff", title=NULL) {
  if (!requireNamespace("ggplot2",quietly=TRUE))
    stop("Requires ggplot2.")
  df <- field_to_df(result,inside_only=TRUE)
  p  <- result$params; L_str <- if(is.finite(p$diffusion_length)) round(p$diffusion_length,3) else "Inf"
  tit <- if(is.null(title)) sprintf("Concentration | %s | D=%g, lambda=%g, L=%s",
                                     p$method,p$D,p$lambda,L_str) else title
  fc  <- if (log_scale){df$lf<-log1p(df$field);"lf"} else "field"
  g   <- ggplot2::ggplot(df,ggplot2::aes(x=x,y=y,fill=.data[[fc]]))+
    ggplot2::geom_raster(interpolate=interpolate)+
    ggplot2::scale_fill_viridis_c(option=palette,
                                   name=if(log_scale)"log1p(C)" else "Conc.",
                                   na.value="transparent")+
    ggplot2::coord_equal()+
    ggplot2::labs(title=tit,x="X",y="Y")+
    ggplot2::theme_minimal(base_size=11)+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold"))
  if (show_contours)
    g <- g+ggplot2::geom_contour(data=df,ggplot2::aes(z=field),
                                  color="white",alpha=0.4,linewidth=0.3,bins=n_contours)
  if (show_sources && !is.null(transcript_coords))
    g <- g+ggplot2::geom_point(data=transcript_coords,ggplot2::aes(x=x,y=y),
                                inherit.aes=FALSE,color=pt_color,size=pt_size,alpha=pt_alpha)
  g
}

sweep_diffusion_length <- function(L_values, mask, transcript_coords, ...) {
  do.call(rbind, lapply(L_values, function(Lv) {
    cat(sprintf("\n--- Sweeping L = %g ---\n",Lv))
    res <- estimate_concentration_field(mask=mask,transcript_coords=transcript_coords,
                                         diffusion_length=Lv,verbose=FALSE,...)
    df  <- field_to_df(res,inside_only=TRUE); df$L <- Lv; df
  }))
}
