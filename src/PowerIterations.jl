
module PowerIterations

	export poweriter
	export orthogonaliter

	using LinearAlgebra

	function poweriter(A::AbstractMatrix{<: T}, iter::Integer) where T <: Number
		poweriter!(A, ones(T, size(A, 1)), iter)
	end
	function poweriter!(A::AbstractMatrix{<: T}, v::AbstractVector{<: T}, iter::Integer) where T <: Number
		@assert ishermitian(A)
		@assert size(A, 1) == size(A, 2) == length(v)
		w = similar(v)
		for i = 1:iter
			normalize!(v)
			mul!(w, A, v)
			v, w = w, v
		end
		λ = Float64(w' * v)
		v ./= norm(v)
		return v, λ
	end

	function orthogonaliter(A::AbstractMatrix{<: T}, rank::Integer, iter::Integer) where T <: Number
		Q = diagm(size(A, 2), rank, ones(T, rank))
		orthogonaliter!(A, Q, rank, iter)
	end
	function orthogonaliter!(A::AbstractMatrix{<: T}, Q::AbstractMatrix{<: T}, rank::Integer, iter::Integer) where T <: Number
		@assert ishermitian(A)
		@assert size(A, 1) == size(A, 2) == size(Q, 1)
		@assert size(A, 2) >= rank
		@assert size(Q, 2) == rank
		V = similar(Q)
		mul!(V, A, Q)
		for i = 1:iter
			F = qr(V)
			mul!(V, A, Matrix(F.Q))
		end
		F = qr(V)
		V = F.Q
		λ = convert.(Float64, diag(F.R))
		V ./= sqrt.(sum(abs2, V; dims=1))
		return V, dropdims(λ; dims=1)
	end

end

