module m_rescaler
    use m_unit_converter
    implicit none

    private
    public t_Rescaler

    !> dx : float
    !>     Grid length [m]
    !> to_c : float
    !>     Light Speed in EMSES
    !> pi : UnitTranslator
    !>     Circular constant []
    !> e : UnitTranslator
    !>     Napiers constant []
    !> c : UnitTranslator
    !>     Light Speed [m/s]
    !> e0 : UnitTranslator
    !>     FS-Permttivity [F/m]
    !> m0 : UnitTranslator
    !>     FS-Permeablity [N/A^2]
    !> qe : UnitTranslator
    !>     Elementary charge [C]
    !> me : UnitTranslator
    !>     Electron mass [kg]
    !> mi : UnitTranslator
    !>     Proton mass [kg]
    !> qe_me : UnitTranslator
    !>     Electron charge-to-mass ratio [C/kg]
    !> kB : UnitTranslator
    !>     Boltzmann constant [J/K]
    !> length : UnitTranslator
    !>     Sim-to-Real length ratio [m]
    !> m : UnitTranslator
    !>     Mass [kg]
    !> t : UnitTranslator
    !>     Time [s]
    !> f : UnitTranslator
    !>     Frequency [Hz]
    !> v : UnitTranslator
    !>     Velocity [m/s]
    !> n : UnitTranslator
    !>     Number density [/m^3]
    !> N : UnitTranslator
    !>     Flux [/m^2s]
    !> F : UnitTranslator
    !>     Force [N]
    !> P : UnitTranslator
    !>     Power [W]
    !> W : UnitTranslator
    !>     Energy [J]
    !> w : UnitTranslator]
    !>     Energy density [J/m^3]
    !> eps : UnitTranslator
    !>     Permittivity  [F/m]
    !> q : UnitTranslator
    !>     Charge [C]
    !> rho : UnitTranslator
    !>     Charge density [C/m^3]
    !> q_m : UnitTranslator
    !>     Charge-to-mass ratio [C/kg]
    !> i : UnitTranslator
    !>     Current [A]
    !> J : UnitTranslator
    !>     Current density [A/m^2]
    !> phi : UnitTranslator
    !>     Potential [V]
    !> E : UnitTranslator
    !>     Electric field [V/m]
    !> H : UnitTranslator
    !>     Magnetic field [A/m]
    !> C : UnitTranslator
    !>     Capacitance [F]
    !> R : UnitTranslator
    !>     Resistance [Ω]
    !> G : UnitTranslator
    !>     Conductance [S]
    !> mu : UnitTranslator
    !>     Permiability [H/m]
    !> B : UnitTranslator
    !>     Magnetic flux density [T]
    !> L : UnitTranslator
    !>     Inductance [H]
    !> T : UnitTranslator
    !>     Temperature [K]
    type :: t_Rescaler
        type(t_UnitConverter) :: dx
        type(t_UnitConverter) :: to_c
        type(t_UnitConverter) :: pi
        type(t_UnitConverter) :: e
        type(t_UnitConverter) :: c
        type(t_UnitConverter) :: e0
        type(t_UnitConverter) :: m0
        type(t_UnitConverter) :: q_m
        type(t_UnitConverter) :: kB
        type(t_UnitConverter) :: length
        type(t_UnitConverter) :: m
        type(t_UnitConverter) :: t
        type(t_UnitConverter) :: f
        type(t_UnitConverter) :: v
        type(t_UnitConverter) :: n
        type(t_UnitConverter) :: flux
        type(t_UnitConverter) :: force
        type(t_UnitConverter) :: P
        type(t_UnitConverter) :: W
        type(t_UnitConverter) :: epsilon
        type(t_UnitConverter) :: q
        type(t_UnitConverter) :: rho
        type(t_UnitConverter) :: i
        type(t_UnitConverter) :: J
        type(t_UnitConverter) :: phi
        type(t_UnitConverter) :: ef
        type(t_UnitConverter) :: bf
        type(t_UnitConverter) :: cap
        type(t_UnitConverter) :: R
        type(t_UnitConverter) :: G
        type(t_UnitConverter) :: mu
        type(t_UnitConverter) :: B
        type(t_UnitConverter) :: L
        type(t_UnitConverter) :: temperature
    end type
end module
