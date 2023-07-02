module m_ohparticles
    use oh_type, only: oh_particle
    implicit none

    private
    public t_OhParticles, new_OhParticles

    type :: t_OhParticles
        type(oh_particle), allocatable :: pbuf(:)
        integer, allocatable :: pbase(:)
        integer(kind=8) :: max_nparticles
        integer :: max_local_particles

        integer, allocatable :: particle_count_histgram(:, :, :) ! (nprocs, nspec, 2)
        integer, allocatable :: total_local_particles(:, :) ! (nspec, 2)

        integer :: nspecies

    contains

        procedure :: allocate_pbuf => ohparticles_allocate_pbuf
    end type

contains

    function new_OhParticles(nspecies, &
                             max_nparticles, &
                             nprocs) result(self)
        integer, intent(in) :: nspecies
        integer(kind=8), intent(in) :: max_nparticles
        integer, intent(in) :: nprocs
        type(t_OhParticles) :: self

        double precision, pointer :: p

        self%nspecies = nspecies
        self%max_nparticles = max_nparticles

        allocate (self%pbase(3))
        allocate (self%particle_count_histgram(nprocs, nspecies, 2))
        allocate (self%total_local_particles(nspecies, 2))
    end function

    subroutine ohparticles_allocate_pbuf(self, max_local_particles)
        class(t_OhParticles), intent(inout) :: self
        integer, intent(in) :: max_local_particles

        self%max_local_particles = max_local_particles
        allocate (self%pbuf(self%max_local_particles))
    end subroutine

    function start_index(self, ispec, ps) result(ret)
        class(t_OhParticles), intent(in) :: self
        integer, intent(in) :: ispec
        integer, intent(in) :: ps
        integer :: ret

        ret = self%pbase(ps) + 1 + sum(self%total_local_particles(1:ispec - 1, ps))
    end function

    function end_index(self, ispec, ps) result(ret)
        class(t_OhParticles), intent(in) :: self
        integer, intent(in) :: ispec
        integer, intent(in) :: ps !> 1: primary subdomain, 2: secondary subdomain
        integer :: ret

        ret = self%pbase(ps) + 1 + sum(self%total_local_particles(1:ispec, ps))
    end function

end module
