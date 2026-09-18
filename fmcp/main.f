      program f_mcp

      use fmcp_data
      use fmcp_read
      use fmcp_interactions
      use fmcp_targets
      use fmcp_flexibility
      use fmcp_dijkstra
      use fmcp_vmd

      implicit none

      call read_trajectories()

      call calculate_interactions()

      call calculate_targets()

      call calculate_flexibility()

      call calculate_dijkstra()

      call write_vmd_tcl()

      end program f_mcp
