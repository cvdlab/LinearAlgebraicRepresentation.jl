function triangulate2D(m::CAGD.Model)
    return Lar.triangulate2D(convert(Lar.Points, m.G'), [m.T[1], m.T[2]])
end