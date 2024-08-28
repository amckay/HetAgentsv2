# # Usage example:
# npts = 10
# rho = 0.9
# sig = 0.01
# gridpts, transition = rouwenhorst(npts, rho, sig)


function RouwenNormalize!(A)
    A[2:end-1, :] ./= 2
    nothing
end

function rouwenupdate!(Πold, Πnew, p, q)
    Nold = size(Πold, 1)
    Πnew[1:Nold, 1:Nold] = p*Πold
    Πnew[1:Nold, 2:(Nold+1)] .+= (1-p)*Πold
    Πnew[2:(Nold+1), 1:Nold] .+= (1-q)*Πold
    Πnew[2:(Nold+1), 2:(Nold+1)] .+= q*Πold
    nothing
end

function rouwenmat(N, p, q)
    initmat = [p 1-p;1-q q]
    if N==2
        ret = initmat
    else
        currmat = initmat
        for n = 3:N
            nextmat = fill(0., (n, n))
            rouwenupdate!(currmat, nextmat, p, q)
            currmat = nextmat
            RouwenNormalize!(currmat)
        end
        ret = currmat
    end
    return(ret)
end

function rouwenhorst(npts, ρ, σ)
    # Discretizes
    # y′ = ρ y +  σ ϵ

    ω = σ/sqrt(1-ρ^2) # Long run std dev
    # println("Long run variance: $ω")
    q = (1+ρ)/2

    points = range(-ω*sqrt(npts - 1),
                   ω*sqrt(npts-1),
                   length=npts)

    # points = GenGrid(uniform(npts, bds))
    Π = rouwenmat(npts, q, q)
    return points, Π
end
