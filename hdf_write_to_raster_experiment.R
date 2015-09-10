library(gdalUtils)
library(raster)
#http://geoscripting-wur.github.io/IntroToRaster/hdf.html



# Get a list of sds names
hdf_file = 'R:/Mann Research/India Night Time Lights/2015 v2/NPP_VDNE_L1.A2015001.2005.P1_03144.2015058164629.hdf'
gdalinfo(hdf_file)

sds <- get_subdatasets('R:/Mann Research/India Night Time Lights/2015 v2/NPP_VDNE_L1.A2015001.2005.P1_03144.2015058164629.hdf')
sds

# Isolate the name of DNB
name <- sds[13]   # radiance 13, lat 15, lon 16
filename <- 'C:/Users/mmann/Desktop/dnb.tif'
gdal_translate(name, dst_dataset = filename)
# Load the Geotiff created into R
dnb = raster(filename)
dnb
plot(dnb)

# Isolate the name of lat
name = sds[15]   # radiance 13, lat 15, lon 16
filename = 'C:/Users/mmann/Desktop/lat.tif'
gdal_translate(name, dst_dataset = filename)
# Load the Geotiff created into R
lat = raster(filename)
lat
plot(lat)

# Isolate the name of lon
name = sds[16]   # radiance 13, lat 15, lon 16
filename = 'C:/Users/mmann/Desktop/lon.tif'
gdal_translate(name, dst_dataset = filename)
# Load the Geotiff created into R
lon = raster(filename)
lon
plot(lon)



a = data.frame(lon=lon[],lat=lat[],dnb=dnb[])
dim(a)

b = rasterFromXYZ(xyz =a,crs='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs',digits=8)
plot(b)




function (xyz, res = c(NA, NA), crs = NA, digits = 5) 
{
    if (length(res) == 1) 
        res = c(res, res)
    if (inherits(xyz, "SpatialPoints")) {
        if (inherits(xyz, "SpatialPointsDataFrame")) {
            xyz <- cbind(coordinates(xyz), xyz@data[, 1])
        }
        else {
            xyz <- coordinates(xyz)
        }
    }
    ln <- colnames(xyz)
    if (inherits(xyz, "data.frame")) {
        xyz <- as.matrix(xyz)
        xyz <- matrix(as.numeric(xyz), ncol = ncol(xyz), nrow = nrow(xyz))
    }
    x <- sort(unique(xyz[, 1]))
    dx <- x[-1] - x[-length(x)]
    if (is.na(res[1])) {
        rx <- min(dx)
        for (i in 1:5) {
            rx <- rx/i
            q <- sum(round(dx/rx, digits = digits)%%1)
            if (q == 0) {
                break
            }
        }
        if (q > 0) {
            stop("x cell sizes are not regular")
        }
    }
    else {
        rx <- res[1]
        test <- sum(round(dx/rx, digits = digits)%%1)
        if (test > 0) {
            stop("x cell sizes are not regular")
        }
    }
    y <- sort(unique(xyz[, 2]))
    dy <- y[-1] - y[-length(y)]
    dy <- round(dy, digits)
    if (is.na(res[2])) {
        ry <- min(dy)
        for (i in 1:5) {
            ry <- ry/i
            q <- sum(round(dy/ry, digits = digits)%%1)
            if (q == 0) {
                break
            }
        }
        if (q > 0) {
            stop("y cell sizes are not regular")
        }
    }
    else {
        ry <- res[2]
        test <- sum(round(dy/ry, digits = digits)%%1)
        if (test > 0) {
            stop("y cell sizes are not regular")
        }
    }
    minx <- min(x) - 0.5 * rx
    maxx <- max(x) + 0.5 * rx
    miny <- min(y) - 0.5 * ry
    maxy <- max(y) + 0.5 * ry
    d <- dim(xyz)
    if (d[2] <= 3) {
        r <- raster(xmn = minx, xmx = maxx, ymn = miny, ymx = maxy, 
                    crs = crs)
    }
    else {
        r <- brick(xmn = minx, xmx = maxx, ymn = miny, ymx = maxy, 
                   crs = crs, nl = d[2] - 2)
    }
    res(r) <- c(rx, ry)
    cells <- cellFromXY(r, xyz[, 1:2])
    if (d[2] > 2) {
        names(r) <- ln[-c(1:2)]
        r[cells] <- xyz[, 3:d[2]]
    }
    return(r)
}






