tags$style(HTML("
  .navbar {
    background: linear-gradient(-45deg, #667eea, #764ba2, #f093fb, #f5576c);
    background-size: 400% 400%;
    animation: gradient 15s ease infinite;
    position: relative;
    overflow: hidden;
  }
  
  @keyframes gradient {
    0% { background-position: 0% 50%; }
    50% { background-position: 100% 50%; }
    100% { background-position: 0% 50%; }
  }
  
  /* Floating particles effect */
  .navbar::before {
    content: '';
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    background-image: 
      radial-gradient(2px 2px at 20px 30px, rgba(255,255,255,0.3), transparent),
      radial-gradient(2px 2px at 40px 70px, rgba(255,255,255,0.3), transparent),
      radial-gradient(1px 1px at 90px 40px, rgba(255,255,255,0.3), transparent),
      radial-gradient(1px 1px at 130px 80px, rgba(255,255,255,0.3), transparent);
    background-repeat: repeat;
    background-size: 200px 100px;
    animation: float 20s linear infinite;
  }
  
  @keyframes float {
    0% { transform: translateX(-200px); }
    100% { transform: translateX(100vw); }
  }
  
  /* Center the title */
  .navbar .container-fluid {
    justify-content: center !important;
  }
  
  .navbar-brand {
    color: white !important;
    font-weight: bold;
    text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
    z-index: 10;
    position: relative;
    margin: 0 !important;
  }
"))