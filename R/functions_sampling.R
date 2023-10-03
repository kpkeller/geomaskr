## Functions For Sampling Points in Circles
## Intended for use with Geomasking of
## participant locations
##
## Author: Kayleigh Keller

#' sample_circle
#'
#' @param x First coordinate for the location(s) to be masked
#' @param y Second coordinate for the location(s) to be mased
#' @param r Radius of masking region
#' @param r0 Minimum radius for masking. Values greater than zero lead to "donut geomasking"
#'
#' @return Data frame containing `x` and `y` coordinates of the perturbed locations
#' @export
#'
#' @examples
#' sample_circle(x=0, y=0, r=1, r0=0)
sample_circle <- function(x, y, r, r0=0){
    n <- length(x)
    r_samp = sqrt(r0^2 + (r^2-r0^2)*(stats::runif(n)))
    theta_samp = stats::runif(n) * 2 * pi
    x = x + r_samp * cos(theta_samp)
    y = y + r_samp * sin(theta_samp)
    data.frame(x=x, y=y)
}



# Version that limits all points to a bounding box,
# and re-samples if outside of it.
# Currently limited to rectangle
#' Title
#'
#' @rdname sample_circle
#' @param xmin lower limit of bounding box
#' @param xmax upper limit of bounding box
#' @param ymin lower limit of bounding box
#' @param ymax upper limit of bounding box
#' @param max_redraw number of times to draw again, if first draw is outside the bounding box
#'
#' @export
#'
sample_circle_boundedbox <- function(x, y, r, r0=0,
                                  xmin=0,
                                  xmax=1,
                                  ymin=0,
                                  ymax=1,
                                  # xrange=c(0, 1),
                                  # yrange=c(0, 1),
                                  max_redraw=100){
    if (length(xmin)==1){
        xmin <- rep(xmin, times=length(x))
    }
    if (length(xmax)==1){
        xmax <- rep(xmax, times=length(x))
    }
    if (length(ymin)==1){
        ymin <- rep(ymin, times=length(x))
    }
    if (length(ymax)==1){
        ymax <- rep(ymax, times=length(x))
    }
    if (length(r)==1){
        r <- rep(r, times=length(x))
    }

    new <- sample_circle(x=x, y=y,
                         r=r, r0=r0)
    inds <- which(new$x < xmin | new$x>xmax |  new$y < ymin | new$y>ymax)
    if (length(inds)>0){
    for (i in 1:length(inds)){
        outside <- TRUE
        for (j in 1:max_redraw){
            # cat("j=", j, "\n")
            new2 <- sample_circle(x=x[inds[i]], y=y[inds[i]],
                                  r=r[inds[i]], r0=r0)
            outside <- any(new2$x < xmin[inds[i]] | new2$x>xmax[inds[i]] |  new2$y < ymin[inds[i]] | new2$y>ymax[inds[i]])
            if(!outside) break();
            if(j==max_redraw) warning("reached max redraw without being in bounds")
        }
        new[inds[i],] <- new2
    }
    }
    new
}





#' @rdname st_sample_radius_bounded
#' @export
#' @importFrom sf st_sample st_coordinates st_crs st_set_crs
#'
#' @details
#' New point is transformed to CRS of `pt`
#'
#'
st_sample_radius <- function(pt, radius=100, region){

    new_pt <- sf::st_sample(region,
                        type="unifdisc",
                        radius=radius,
                        centre=sf::st_coordinates(pt),
                        n=1)
    new_pt <- sf::st_set_crs(new_pt, sf::st_crs(pt))
    new_pt
}

#' Title
#' Sample Points within a Radius
#'
#' @description Randomly samples points within a given radius using  functions from `sf` and `spatstat`.
#'
#' @param pt Point location to perturb. Assumed to be an
#' @param radius Radius of circle to sample within
#' @param region `sfc` object
#' @param minradius Minimum radius
#' @param maxretry Number of times to retry point selection if selected point is either inside minimum radius or not in bounding region
#' @param return_dist Logical indicator of whether or not to return the distance of new point from original point
#'
#' @return an `sfc` object containg the new point
#' @export
#' @importFrom sf st_intersects st_distance
#'
# #' @examples
st_sample_radius_bounded <- function(pt,
                                       radius=100,
                                       region,
                                       minradius=0,
                                       maxretry=1000,
                                       return_dist=TRUE){

    new_pt <- st_sample_radius(pt,
                                 radius=radius,
                                 region=region)
    inside <- sf::st_intersects(region, new_pt, sparse=F)
    move_dist <-  as.numeric(sf::st_distance(pt, new_pt))
    if (minradius>0){
        tooclose <- move_dist < minradius
    } else {
        tooclose <- FALSE
    }


    if((!inside) | tooclose){
        for (i in 1:maxretry){
            newer_pt <- st_sample_radius(pt,
                                           radius=radius,
                                           region=region)
            move_dist <-  as.numeric(st_distance(pt, newer_pt))
            if(st_intersects(region, newer_pt, sparse=F) & move_dist > minradius) break;
        }
        if (i==maxretry) warning("max resample iteration reached.")
        new_pt <- newer_pt
    }
    if (!return_dist){
        out <- new_pt
    } else {
        out <-  list(pt=new_pt,
                  dist=move_dist)
    }
    out
}

##' @rdname st_sample_radius_bounded
##' @param pts Collection of points to resample location for
##' @export
st_sample_radius_bounded_set <- function(pts,
                                           radius=100,
                                           region,
                                           minradius=0,
                                           maxretry=1000){
    newpts <- data.frame()
    newdists <- numeric()
    for (j in 1:length(pts)){
        newpts_temp <- st_sample_radius_bounded(pt=pts[j],
                                                  radius=radius,
                                                  region=region,
                                                  minradius=minradius,
                                                  maxretry=maxretry,
                                                  return_dist=TRUE)
        newpts <- rbind(newpts,
                        newpts_temp$pt)
        newdists <- rbind(newdists,
                          newpts_temp$dist)
    }
    list(pt=newpts,
         dist=newdists)
}
