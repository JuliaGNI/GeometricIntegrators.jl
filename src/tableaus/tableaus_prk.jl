
function CoefficientsSymplecticEulerForward()
    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [1.0, 0.0]
    c = [0.0, 1.0]

    o = 1

    Tableau(:symplectic_euler_forward, o, a, b, c)
end

function CoefficientsSymplecticEulerBackward()
    a = [[0.0 0.0]
         [0.0 1.0]]
    b = [0.0, 1.0]
    c = [0.0, 1.0]

    o = 1

    Tableau(:symplectic_euler_backward, o, a, b, c)
end

"Tableau for symplectic Euler-A method"
function TableauSymplecticEulerA()
    TableauPRK(:symplectic_euler_a, CoefficientsSymplecticEulerForward(), CoefficientsSymplecticEulerBackward())
end

"Tableau for symplectic Euler-B method"
function TableauSymplecticEulerB()
    TableauPRK(:symplectic_euler_b, CoefficientsSymplecticEulerBackward(), CoefficientsSymplecticEulerForward())
end


"Gauss-Legendre Runge-Kutta"
function PartitionedTableauGauss(s::Int)
    TableauPRK(Symbol("PGauss", s), TableauGauss(s))
end


"Tableau for Gauss-Lobatto IIIAIIIB method with s stages"
function TableauLobattoIIIAIIIB(s)
    TableauPRK(:lobatto_IIIA_IIIB, TableauLobattoIIIA(s), TableauLobattoIIIB(s))
end

"Tableau for Gauss-Lobatto IIIBIIIA method with s=2 stages"
function TableauLobattoIIIBIIIA(s)
    TableauPRK(:lobatto_IIIB_IIIA, TableauLobattoIIIB(s), TableauLobattoIIIA(s))
end
