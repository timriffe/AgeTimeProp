# help functions for diagram of lifeline
# clockwise quarter arc (90 degrees)
degrees2radians <- function(degrees){
    degrees * (pi / 180)
}

quarterArc <- function(x, y, radius = 1, fromDegrees = 180, ...){
    xx <- degrees2radians(seq(fromDegrees, fromDegrees + 90, by = .5))
    x <- cos(xx) * radius + x
    y <- sin(xx) * radius + y
    lines(x, y, ...)
}

curlyBrace1 <- function(xl, y, length = 5, radius1 = .5, radius2 = .25, top = TRUE, ...){  
    # i.e. the pointy part on top or on bottom?
    if (top){
        quarterArc(xl + radius1, y - radius1, radius = radius1, fromDegrees = 90, ...)
        quarterArc(xl + length - radius1, y - radius1 , radius = radius1, fromDegrees = 0, ...)
        quarterArc(xl + length / 2 - radius2, y + radius2, radius = radius2, fromDegrees = 270, ...)
        quarterArc(xl + length / 2 + radius2, y + radius2, radius = radius2, fromDegrees = 180, ...)
    } else {
        quarterArc(xl + radius1, y + radius1, radius = radius1, fromDegrees = 180, ...)
        quarterArc(xl + length - radius1, y + radius1 , radius = radius1, fromDegrees = 0 - 90, ...)
        quarterArc(xl + length / 2 - radius2, y - radius2, radius = radius2, fromDegrees = 270 + 90, ...)
        quarterArc(xl + length / 2 + radius2, y - radius2, radius = radius2, fromDegrees = 180 - 90, ...)        
    }
    segments(xl + radius1, y, xl + length / 2 - radius2, y, ...)
    segments(xl + length - radius1, y, xl + length / 2 + radius2, y, ...)   
}
