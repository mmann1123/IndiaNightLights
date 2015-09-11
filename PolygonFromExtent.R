PolygonFromExtent <-
function(ext, asSpatial=T, crs=CRS(NA), id=1)
{
  if(class(ext)== "rasterLayer")
  {
    # if raster supplied determine extent and crs then proceed
    crs <- ext@crs
    ext <- extent(ext)
  }
  if(class(ext)== "Extent")
        x1 <- ext@xmin
        x2 <- ext@xmax
        y1 <- ext@ymin
        y2<-ext@ymax

        coords <- matrix(c(x1, y1,
                                           x1, y2,
                                           x2, y2,
                                           x2, y1,
                                           x1, y1), ncol=2, byrow=T)

        poly <- Polygon(coords)
        if(asSpatial)
        {
                spPoly <- SpatialPolygons(list(Polygons(list(poly), ID=id)), proj4string=crs)
                return(spPoly)

        }
        return(poly)

}

