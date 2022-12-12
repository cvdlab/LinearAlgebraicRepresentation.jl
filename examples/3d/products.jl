using LinearAlgebraicRepresentation
using Plasm
using DataStructures
Lar = LinearAlgebraicRepresentation

v0 = [0. 1. 2. 3. 4. 5. 6. 7. 8.]
v1 = collect(0.:12)
c0 = Lar.larGridSkeleton([1,1,1])(0)
m0 = (v0, c0)
c1 = [[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]]
v1 = hcat(collect(0.:12)...)
m0 = (v0,c0)
m1 = (v1,c1)
Plasm.view(Lar.larModelProduct( m1,m0 ))

m11 = Lar.larModelProduct( m0, m1 )
m12 = Lar.larModelProduct( m1, m0 )
Lar.Struct([ m11,m12 ])

m2 = Lar.struct2lar(Lar.Struct([ m11,m12 ]))
Plasm.view(m2)

m21 = Lar.larModelProduct( m2,m1 )
Plasm.view(m21)

m20 = Lar.larModelProduct( m1,m1 )
Plasm.view(m20)

m20 = Lar.larModelProduct( m20,m0 )
m20 = Lar.larModelProduct( m1,m1 )
m22 = Lar.larModelProduct( m20,m0 )
Plasm.view(m22)

Lar.struct2lar(Lar.Struct([ m22,m21 ]))
m32 = Lar.struct2lar(Lar.Struct([ m22,m21 ]))
Plasm.view(m32)










