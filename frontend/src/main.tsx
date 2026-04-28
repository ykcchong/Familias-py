import React from "react";
import ReactDOM from "react-dom/client";
import { BrowserRouter, Routes, Route } from "react-router-dom";
import App from "./App";
import HomePage from "./pages/HomePage";
import OneParentPage from "./pages/OneParentPage";
import TwoParentPage from "./pages/TwoParentPage";
import ArbitraryPage from "./pages/ArbitraryPage";
import "ag-grid-community/styles/ag-grid.css";
import "ag-grid-community/styles/ag-theme-quartz.css";
import "./styles.css";

ReactDOM.createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    <BrowserRouter>
      <Routes>
        <Route element={<App />}>
          <Route index element={<HomePage />} />
          <Route path="one-parent" element={<OneParentPage />} />
          <Route path="two-parent" element={<TwoParentPage />} />
          <Route path="arbitrary" element={<ArbitraryPage />} />
        </Route>
      </Routes>
    </BrowserRouter>
  </React.StrictMode>
);
