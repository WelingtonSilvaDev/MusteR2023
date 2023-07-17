#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= multiple viewer =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#mol3d_m("7q9y") #exemplo de chamada
mol3d_multi <- function(pdb_view_m){
  library(r3dmol)
  # Set up the initial viewer
  r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 10,
      lowerZoomLimit = 50,
      upperZoomLimit = 350
    ),
  ) %>%
  m1 <- r3dmol() %>%
    m_add_model(data = paste0("pdb/", pdb_view_m, ".pdb"), format = "pdb") %>%
    m_add_model(data = pdb_1j72, format = "pdb") %>%
    m_zoom_to()
  
  m2 <- m1 %>%
    m_set_style(style = m_style_cartoon(color = "spectrum"))
  
  m3 <- r3dmol() %>%
    m_add_model(data = paste0("pdb/", pdb_view_m, ".pdb"), format = "pdb") %>%
    m_set_style(style = m_style_stick()) %>%
    m_zoom_to()
  
  m4 <- m3 %>%
    m_set_style(style = m_style_sphere())
  
  m_grid(
    viewer = list(m1, m2, m3, m4),
    rows = 2,
    cols = 2,
    control_all = TRUE,
    viewer_config = m_viewer_spec(
      backgroundColor = "lightblue"
    )
  )
}
