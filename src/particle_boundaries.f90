module m_particle_boundaries
    use m_particle_boundary
    use oh_type, only: oh_particle
    use m_periodic_particle_boundary
    implicit none

    private
    public t_ParticleBoundaries
    public new_ParticleBoundaries

    type :: t_ParticleBoundaryHolder
        class(t_ParticleBoundary), allocatable :: particle_boundary
    end type

    type, extends(t_ParticleBoundary) :: t_ParticleBoundaries
        type(t_ParticleBoundaryHolder) :: particle_boundaries(3)
    contains
        procedure :: apply => particleBoundaries_apply
    end type

contains

    function new_ParticleBoundaries(nx, ny, nz) result(obj)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        type(t_ParticleBoundaries) :: obj

        obj%particle_boundaries(1)%particle_boundary = new_PeriodicParticleBoundary(1, nx)
        obj%particle_boundaries(2)%particle_boundary = new_PeriodicParticleBoundary(2, ny)
        obj%particle_boundaries(3)%particle_boundary = new_PeriodicParticleBoundary(3, nz)
    end function

    subroutine particleBoundaries_apply(self, particle, dt)
        class(t_ParticleBoundaries), intent(in) :: self
        type(oh_particle), intent(inout) :: particle
        double precision, intent(in) :: dt
    
        integer :: i
        do i = 1, size(self%particle_boundaries)
            call self%particle_boundaries(i)%particle_boundary%apply(particle, dt)
        end do
    end subroutine
    
end module
