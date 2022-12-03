using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

# using PyCall
# Io3D = @pyimport("io3d")

function randomImage(shape, structure, noiseFraction=0.0)
   """ Generation of random image of given shape and structure. 
      Return scipy.ndarray(shape)
   """
   println("noiseFraction = $noiseFraction")
   ranges = [shape[k]/structure[k] for k=1:length(shape)]
   println("ranges = $ranges")
   random_array = rand(Int8, structure)
end

shape = [4,4]
structure = [8,8]
random_array = randomImage(shape, structure, 0.0)


function boundaryOps(skeletons)

	characteristicMatrices = map(Lar.characteristicMatrix,skeletons)
	couples = zip(characteristicMatrices[1:end-1], characteristicMatrices[2:end])
	operators = [ div(,2) for (A,B) in couples ]

end


function imageChainComplex (shape)
	_,skeletons = Lar.cuboidGrid(shape,true)
	
	
	
	
		
	
	imageLAR = (shape, skeletons, operators)
	dump(imageLAR,filename)
	print "filename =",filename
	return shape, skeletons, operators
end



def imageChainComplex (shape):
   tokens = str(shape)[1:-1].split(',')
   tokens = [token.strip() for token in tokens]
   filename = "tmp/larimage-" + "-".join(tokens) + ".pickle"
   
   if os.path.isfile(filename):
      shape, skeletons, operators = loadImageLAR(filename)
   else:
      skeletons = gridSkeletons(list(shape))
      operators = boundaryOps(skeletons)
      imageLAR = (shape, skeletons, operators)
      dump(imageLAR,filename)
      print "filename =",filename
   return shape, skeletons, operators
